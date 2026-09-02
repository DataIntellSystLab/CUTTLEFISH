import java.io.*;
import java.math.BigInteger;
import java.lang.reflect.Method;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import com.clearspring.analytics.hash.MurmurHash;

public class ShardedBloomFilter
{
    private static final ThreadLocal<byte[]> LONG_BUFFER =
        ThreadLocal.withInitial(() -> new byte[8]);

    private BloomFilterSuite[] filters;

    // Single-writer reusable workspace for shard-grouped population batches.
    // Construction/population is sequential per filter, so no synchronization is needed.
    private long[] groupedHashWorkspace = new long[0];
    private int[] shardCountsWorkspace;
    private int[] shardOffsetsWorkspace;
    // Reused counting-sort cursor array; avoids Arrays.copyOf per super-batch.
    private int[] shardNextWorkspace;
    private int numShards;
    private String baseFilename;

    public ShardedBloomFilter(long maxElements, double falsePositiveRate, String filename)
    {
        this.baseFilename = filename;

        long maxCapacityPerShard = getMaxCapacityPerShard(falsePositiveRate);
        int requiredShards = (int)Math.ceil(
            (double)maxElements / (double)maxCapacityPerShard
        );

        numShards = 1;
        while (numShards < requiredShards)
            numShards <<= 1;

        if (numShards > 255)
            throw new IllegalArgumentException(
                "The number of shards cannot exceed 255 because it is stored in one header byte"
            );

        long capacityPerShard = (long)Math.ceil(
            (double)maxElements / (double)numShards
        );

        if (capacityPerShard > Integer.MAX_VALUE)
            throw new IllegalArgumentException(
                "Capacity per shard exceeds the supported integer range"
            );

        filters = new BloomFilterSuite[numShards];
        for (int i = 0; i < numShards; i++)
        {
            filters[i] = new BloomFilterSuite(
                (int)capacityPerShard,
                falsePositiveRate,
                filename + "." + i,
                numShards,
                numShards == 1
            );
        }
    }

    public ShardedBloomFilter(String filename, long inMemoryLimitBytes)
    throws Exception
    {
        this.baseFilename = filename;
        this.numShards = BloomFilterSuite.readNumShards(filename + ".0");
        this.filters = new BloomFilterSuite[numShards];

        long totalSize = 0;
        for (int i = 0; i < numShards; i++)
        {
            File f = new File(filename + "." + i);
            if (!f.exists())
                throw new RuntimeException(
                    "Missing shard file: " + f.getAbsolutePath()
                );
            totalSize += f.length();
        }

        boolean loadIntoMemory = totalSize <= Math.max(0L,inMemoryLimitBytes);
        int openMode = loadIntoMemory ? Integer.MAX_VALUE : 0;

        for (int i = 0; i < numShards; i++)
        {
            filters[i] = new BloomFilterSuite();
            filters[i].openReadOnly(filename + "." + i, openMode);

            if (filters[i].getNumShards() != numShards)
                throw new IOException(
                    "Shard-count mismatch in file " + filename + "." + i
                );
        }
    }

    // Returns the seed-0 hash in the high 32 bits and the seed-1 hash in
    // the low 32 bits. The two seed-independent word-mixing operations are
    // performed once, while the Murmur accumulator and final avalanche remain
    // independent for each seed. This is bit-for-bit compatible with hashing
    // the original eight-byte big-endian representation with seeds 0 and 1.
    private static long hashLongBigEndianBothSeeds(long value)
    {
        final int m = 0x5bd1e995;
        final int r = 24;

        int highWord = Integer.reverseBytes((int)(value >>> 32));
        highWord *= m;
        highWord ^= highWord >>> r;
        highWord *= m;

        int lowWord = Integer.reverseBytes((int)value);
        lowWord *= m;
        lowWord ^= lowWord >>> r;
        lowWord *= m;

        int hash0 = Long.BYTES;
        int hash1 = 1 ^ Long.BYTES;

        hash0 *= m;
        hash0 ^= highWord;
        hash1 *= m;
        hash1 ^= highWord;

        hash0 *= m;
        hash0 ^= lowWord;
        hash1 *= m;
        hash1 ^= lowWord;

        hash0 ^= hash0 >>> 13;
        hash0 *= m;
        hash0 ^= hash0 >>> 15;

        hash1 ^= hash1 >>> 13;
        hash1 *= m;
        hash1 ^= hash1 >>> 15;

        return ((long)hash0 << 32) | Integer.toUnsignedLong(hash1);
    }
    private int getShard(byte[] key)
    {
        return MurmurHash.hash(key, 0) & (numShards - 1);
    }

    private static byte[] longToBytes(long value)
    {
        byte[] bytes = LONG_BUFFER.get();
        bytes[0] = (byte)(value >>> 56);
        bytes[1] = (byte)(value >>> 48);
        bytes[2] = (byte)(value >>> 40);
        bytes[3] = (byte)(value >>> 32);
        bytes[4] = (byte)(value >>> 24);
        bytes[5] = (byte)(value >>> 16);
        bytes[6] = (byte)(value >>> 8);
        bytes[7] = (byte)value;
        return bytes;
    }

    public void add(byte[] key)
    {
        if (numShards == 1)
        {
            filters[0].add(key);
            return;
        }
        filters[getShard(key)].add(key);
    }

    public boolean contains(byte[] key) throws Exception
    {
        if (numShards == 1)
            return filters[0].contains(key);
        return filters[getShard(key)].contains(key);
    }

    public void add(long value)
    {
        byte[] bytes = longToBytes(value);
        if (numShards == 1)
        {
            filters[0].add(bytes);
            return;
        }
        filters[getShard(bytes)].add(bytes);
    }
    public void addSingleWriter(long value)
    {
        long hashes = hashLongBigEndianBothSeeds(value);
        int hash = (int)(hashes >>> 32);
        int stepHash = (int)hashes;
        int shard = numShards == 1 ? 0 : hash & (numShards - 1);
        filters[shard].addHashed(
            Integer.toUnsignedLong(hash),
            Integer.toUnsignedLong(stepHash)
        );
    }

    /**
     * Exact single-writer batch insertion with one hash calculation per value.
     * The input batch is converted in place from minimizers to packed hash pairs,
     * then stably grouped by target shard in reusable workspace. Bloom positions
     * and resulting filter bits are identical to addSingleWriter().
     */
    public void addBatchGroupedByShard(long[] values, int length)
    {
        if (length < 0 || length > values.length)
            throw new IllegalArgumentException("Invalid batch length: " + length);
        if (length == 0)
            return;

        if (numShards == 1)
        {
            BloomFilterSuite filter = filters[0];
            for (int i = 0; i < length; i++)
            {
                long hashes = hashLongBigEndianBothSeeds(values[i]);
                filter.addHashed(
                    Integer.toUnsignedLong((int)(hashes >>> 32)),
                    Integer.toUnsignedLong((int)hashes)
                );
            }
            return;
        }

        if (groupedHashWorkspace.length < length)
            groupedHashWorkspace = new long[length];
        if (shardCountsWorkspace == null || shardCountsWorkspace.length != numShards)
        {
            shardCountsWorkspace = new int[numShards];
            shardOffsetsWorkspace = new int[numShards];
            shardNextWorkspace = new int[numShards];
        }
        Arrays.fill(shardCountsWorkspace, 0);

        // Hash once, retain the packed pair in the caller's reusable batch,
        // and count its destination shard.
        for (int i = 0; i < length; i++)
        {
            long hashes = hashLongBigEndianBothSeeds(values[i]);
            values[i] = hashes;
            int shard = ((int)(hashes >>> 32)) & (numShards - 1);
            shardCountsWorkspace[shard]++;
        }

        int position = 0;
        for (int shard = 0; shard < numShards; shard++)
        {
            shardOffsetsWorkspace[shard] = position;
            position += shardCountsWorkspace[shard];
        }

        // Stable counting-sort partition by shard.
        for (int i = 0; i < length; i++)
        {
            long hashes = values[i];
            int shard = ((int)(hashes >>> 32)) & (numShards - 1);
            groupedHashWorkspace[shardOffsetsWorkspace[shard]++] = hashes;
        }

        position = 0;
        for (int shard = 0; shard < numShards; shard++)
        {
            int end = position + shardCountsWorkspace[shard];
            BloomFilterSuite filter = filters[shard];
            for (int i = position; i < end; i++)
            {
                long hashes = groupedHashWorkspace[i];
                filter.addHashed(
                    Integer.toUnsignedLong((int)(hashes >>> 32)),
                    Integer.toUnsignedLong((int)hashes)
                );
            }
            position = end;
        }
    }

    public boolean contains(long value) throws Exception
    {
        byte[] bytes = longToBytes(value);
        if (numShards == 1)
            return filters[0].contains(bytes);
        return filters[getShard(bytes)].contains(bytes);
    }

    public long getNumElements()
    {
        long total = 0;
        for (BloomFilterSuite f : filters)
            total += f.getNumElements();
        return total;
    }

    public long getEstimatedDistinctElements()
    {
        long total = 0;
        for (BloomFilterSuite f : filters)
            total += f.getEstimatedDistinctElements();
        return total;
    }

    public double calculateFalsePositiveRate()
    {
        if(filters.length==0)return 0.0;
        double result=0.0, weight=1.0/filters.length;
        for(BloomFilterSuite filter:filters)result+=weight*filter.calculateFalsePositiveRate();
        return result;
    }

    /**
     * Exact mapped-parallel insertion for multi-shard filters. Values are hashed
     * once, counting-sorted by shard in reusable workspace, and each shard range
     * is updated by one exclusive task. The call returns only after all tasks
     * complete, so the workspace can be reused safely by the next super-batch.
     */
    public void addBatchMappedParallel(long[] values, int length,
            ExecutorService executor) throws Exception
    {
        if (executor == null || numShards == 1)
        {
            addBatchGroupedByShard(values,length);
            return;
        }
        if (length < 0 || length > values.length)
            throw new IllegalArgumentException("Invalid batch length: " + length);
        if (length == 0) return;

        if (groupedHashWorkspace.length < length)
            groupedHashWorkspace = new long[length];
        if (shardCountsWorkspace == null || shardCountsWorkspace.length != numShards)
        {
            shardCountsWorkspace = new int[numShards];
            shardOffsetsWorkspace = new int[numShards];
            shardNextWorkspace = new int[numShards];
        }
        Arrays.fill(shardCountsWorkspace,0);
        int shardMask=numShards-1;
        for (int i=0;i<length;i++)
        {
            long hashes=hashLongBigEndianBothSeeds(values[i]);
            values[i]=hashes;
            shardCountsWorkspace[((int)(hashes>>>32))&shardMask]++;
        }
        int position=0;
        for (int shard=0;shard<numShards;shard++)
        {
            shardOffsetsWorkspace[shard]=position;
            position+=shardCountsWorkspace[shard];
        }
        System.arraycopy(shardOffsetsWorkspace,0,shardNextWorkspace,0,numShards);
        for (int i=0;i<length;i++)
        {
            long hashes=values[i];
            int shard=((int)(hashes>>>32))&shardMask;
            groupedHashWorkspace[shardNextWorkspace[shard]++]=hashes;
        }

        ArrayList<Callable<Void>> tasks=new ArrayList<Callable<Void>>(numShards);
        for (int shard=0;shard<numShards;shard++)
        {
            final int ownedShard=shard;
            final int from=shardOffsetsWorkspace[shard];
            final int to=from+shardCountsWorkspace[shard];
            if (from==to) continue;
            tasks.add(() ->
            {
                BloomFilterSuite filter=filters[ownedShard];
                for (int i=from;i<to;i++)
                {
                    long hashes=groupedHashWorkspace[i];
                    filter.addHashed(Integer.toUnsignedLong((int)(hashes>>>32)),
                        Integer.toUnsignedLong((int)hashes));
                }
                return null;
            });
        }
        for (Future<Void> future:executor.invokeAll(tasks)) future.get();
    }

    public int getNumShards()
    {
        return numShards;
    }

    public long getDiskSizeBytes()
    {
        long total = 0;
        for (int i = 0; i < numShards; i++)
        {
            File f = new File(baseFilename + "." + i);
            if (f.exists())
                total += f.length();
        }
        return total;
    }

    public static long getDiskSizeBytesForBase(String filename) throws IOException
    {
        int shards=BloomFilterSuite.readNumShards(filename+".0");
        long total=0L;
        for(int i=0;i<shards;i++)
        {
            File file=new File(filename+"."+i);
            if(!file.isFile())throw new IOException("Missing shard file: "+file.getAbsolutePath());
            total=Math.addExact(total,file.length());
        }
        return total;
    }

    public double getDiskSizeMB()
    {
        return getDiskSizeBytes() / 1024.0 / 1024.0;
    }

    public void close() throws Exception
    {
        for (BloomFilterSuite f : filters)
            f.close();
    }

    private static long getMaxCapacityPerShard(double falsePositiveRate)
    {
        return (long)(
            Integer.MAX_VALUE * Math.pow(Math.log(2), 2)
            / -Math.log(falsePositiveRate)
        );
    }

    private static class BloomFilterSuite
    {
        private int numBits;
        private int numHashes;
        private long numElements;
        private int numShards;
        public RandomAccessFile file;
        public MappedByteBuffer buffer;
        public byte[] inMemoryFilter;
        private boolean readOnly;
        private String heapBuildFilename;

        // Header layout, from the end of the file:
        // [numElements: 8 bytes][numShards: 1 byte]
        // [remainingBits: 1 byte][numHashes: 1 byte]
        private static final int HEADER_BYTES = 11;
        private static final int NUM_ELEMENTS_OFFSET_FROM_END = 11;
        private static final int NUM_SHARDS_OFFSET_FROM_END = 3;
        private static final int REMAINING_BITS_OFFSET_FROM_END = 2;
        private static final int NUM_HASHES_OFFSET_FROM_END = 1;

        private static final int[] BYTE_COUNTS = new int[256];

        static
        {
            for (int i = 0; i < 256; i++)
                BYTE_COUNTS[i] = Integer.bitCount(i);
        }

        public BloomFilterSuite()
        {
            numBits = 0;
            numHashes = 0;
            numElements = 0;
            numShards = 0;
            file = null;
            buffer = null;
            inMemoryFilter = null;
            readOnly = false;
            heapBuildFilename = null;
        }

        public BloomFilterSuite(
            int capacity,
            double falsePositiveRate,
            String filename,
            int numShards,
            boolean heapBuild
        )
        {
            if (numShards < 1 || numShards > 255)
                throw new IllegalArgumentException(
                    "numShards must be between 1 and 255"
                );

            this.numBits = calculateNumBits(capacity, falsePositiveRate);
            this.numHashes = calculateNumHashes(numBits, capacity);
            this.numElements = 0;
            this.numShards = numShards;
            this.readOnly = false;

            if (numHashes > 255)
                throw new IllegalArgumentException(
                    "Cannot initialize Bloom filter with more than 255 hashes"
                );

            int dataBytes = (int)(((long)numBits + 7L) >>> 3);
            int fileLength = dataBytes + HEADER_BYTES;
            byte[] emptyBytes = new byte[fileLength];
            Arrays.fill(emptyBytes, (byte)0);

            ByteBuffer.wrap(
                emptyBytes,
                fileLength - NUM_ELEMENTS_OFFSET_FROM_END,
                Long.BYTES
            ).putLong(0L);
            emptyBytes[fileLength - NUM_SHARDS_OFFSET_FROM_END] =
                (byte)numShards;
            emptyBytes[fileLength - REMAINING_BITS_OFFSET_FROM_END] =
                (byte)(numBits % 8);
            emptyBytes[fileLength - NUM_HASHES_OFFSET_FROM_END] =
                (byte)(numHashes - 128);

            if (filename.isEmpty() || heapBuild)
            {
                inMemoryFilter = emptyBytes;
                buffer = null;
                file = null;
                heapBuildFilename = filename.isEmpty() ? null : filename;
            }
            else
            {
                try
                {
                    File f = new File(filename);
                    if (f.exists() && !f.delete())
                        throw new IOException(
                            "Could not replace existing file: " + filename
                        );

                    file = new RandomAccessFile(filename, "rw");
                    file.write(emptyBytes);
                    buffer = file.getChannel().map(
                        FileChannel.MapMode.READ_WRITE,
                        0,
                        fileLength
                    );
                    inMemoryFilter = null;
                    heapBuildFilename = null;
                }
                catch (IOException e)
                {
                    throw new RuntimeException(
                        "Error creating memory-mapped file: " + e.getMessage(),
                        e
                    );
                }
            }
        }

        public int getNumBits()
        {
            return numBits;
        }

        public int getNumHashes()
        {
            return numHashes;
        }

        public int getNumShards()
        {
            return numShards;
        }

        public boolean isInMemory()
        {
            return inMemoryFilter != null;
        }

        public static int readNumShards(String filename) throws IOException
        {
            try (RandomAccessFile shard = new RandomAccessFile(filename, "r"))
            {
                long fileLength = shard.length();
                if (fileLength < HEADER_BYTES)
                    throw new IOException(
                        "Invalid Bloom filter shard: " + filename
                    );

                shard.seek(fileLength - NUM_SHARDS_OFFSET_FROM_END);
                int storedNumShards = shard.readUnsignedByte();
                if (storedNumShards < 1)
                    throw new IOException(
                        "Invalid shard count in file: " + filename
                    );
                return storedNumShards;
            }
        }

        public void openReadOnly(String filename, int inMemoryLimitBytes)
        throws Exception
        {
            this.readOnly = true;
            this.file = new RandomAccessFile(filename, "r");
            int fileLength = (int)file.length();

            if (fileLength < HEADER_BYTES)
            {
                file.close();
                file = null;
                throw new IOException("Invalid Bloom filter shard: " + filename);
            }

            this.buffer = file.getChannel().map(
                FileChannel.MapMode.READ_ONLY,
                0,
                fileLength
            );
            int dataBytes = fileLength - HEADER_BYTES;
            int usedBitsInLastByte = Byte.toUnsignedInt(
                buffer.get(fileLength - REMAINING_BITS_OFFSET_FROM_END)
            );
            this.numBits = usedBitsInLastByte == 0
                ? dataBytes * 8
                : (dataBytes - 1) * 8 + usedBitsInLastByte;
            this.numElements = buffer.getLong(
                fileLength - NUM_ELEMENTS_OFFSET_FROM_END
            );
            this.numShards = Byte.toUnsignedInt(
                buffer.get(fileLength - NUM_SHARDS_OFFSET_FROM_END)
            );

            if (numShards < 1)
                throw new IOException(
                    "Invalid shard count in file: " + filename
                );

            if ((long)inMemoryLimitBytes >= fileLength)
            {
                this.inMemoryFilter = new byte[fileLength];
                buffer.position(0);
                buffer.get(this.inMemoryFilter);
                this.numHashes =
                    (int)this.inMemoryFilter[
                        fileLength - NUM_HASHES_OFFSET_FROM_END
                    ] + 128;
                buffer.clear();
                closeDirectBuffer(this.buffer);
                this.file.close();
                this.buffer = null;
                this.file = null;
            }
            else
            {
                this.numHashes =
                    (int)buffer.get(
                        fileLength - NUM_HASHES_OFFSET_FROM_END
                    ) + 128;
                this.inMemoryFilter = null;
            }
        }

        private int calculateNumBits(int capacity, double falsePositiveRate)
        {
            long requested=(long)Math.ceil(
                -(capacity*Math.log(falsePositiveRate))
                /Math.pow(Math.log(2),2)
            );
            if(requested<2L)requested=2L;
            if(requested>Integer.MAX_VALUE)
                throw new IllegalArgumentException("Bloom shard requires more than Integer.MAX_VALUE bits");
            int candidate=(int)requested;
            if(BigInteger.valueOf(candidate).isProbablePrime(20))return candidate;
            BigInteger prime=BigInteger.valueOf(candidate).nextProbablePrime();
            if(prime.compareTo(BigInteger.valueOf(Integer.MAX_VALUE))>0)
                throw new IllegalArgumentException("Unable to size Bloom shard with a 32-bit prime bit count");
            return prime.intValue();
        }

        private int calculateNumHashes(int numBits, int capacity)
        {
            return (int)Math.ceil(
                (numBits / (double)capacity) * Math.log(2)
            );
        }

        public double calculateFalsePositiveRate()
        {
            return Math.pow(
                1.0 - Math.exp(
                    -(double)numHashes * (double)numElements
                    / (double)numBits
                ),
                numHashes
            );
        }

        public synchronized void add(byte[] data)
        {
            long hash = Integer.toUnsignedLong(MurmurHash.hash(data, 0));
            long step = Integer.toUnsignedLong(MurmurHash.hash(data, 1));
            long bitIndex = hash % numBits;
            step %= numBits;
            if (step == 0)
                step = 1; // numBits is prime: every nonzero step spans the full cycle

            for (int i = 0; i < numHashes; i++)
            {
                int bit = (int)bitIndex;
                int pos = bit >> 3;
                int cbit = 1 << (bit & 0x7);
                byte sbit = inMemoryFilter == null
                    ? buffer.get(pos)
                    : inMemoryFilter[pos];
                sbit = (byte)(sbit | cbit);

                if (inMemoryFilter == null)
                    buffer.put(pos, sbit);
                else
                    inMemoryFilter[pos] = sbit;

                bitIndex += step;
                if (bitIndex >= numBits)
                    bitIndex -= numBits;
            }
            numElements++;
        }

        public void addHashed(long hash, long step)
        {
            long bitIndex = hash % numBits;
            step %= numBits;
            if (step == 0)
                step = 1; // numBits is prime: every nonzero step spans the full cycle
            MappedByteBuffer localBuffer = buffer;
            byte[] localFilter = inMemoryFilter;
            int localNumBits = numBits;
            int localNumHashes = numHashes;
            for (int i = 0; i < localNumHashes; i++)
            {
                int bit = (int)bitIndex;
                int pos = bit >>> 3;
                int cbit = 1 << (bit & 0x7);
                if (localFilter == null)
                    localBuffer.put(pos, (byte)(localBuffer.get(pos) | cbit));
                else
                    localFilter[pos] = (byte)(localFilter[pos] | cbit);
                bitIndex += step;
                if (bitIndex >= localNumBits)
                    bitIndex -= localNumBits;
            }
            numElements++;
        }
        public boolean contains(byte[] data) throws Exception
        {
            long hash = Integer.toUnsignedLong(MurmurHash.hash(data, 0));
            long step = Integer.toUnsignedLong(MurmurHash.hash(data, 1));
            long bitIndex = hash % numBits;
            step %= numBits;
            if (step == 0)
                step = 1; // numBits is prime: every nonzero step spans the full cycle

            for (int i = 0; i < numHashes; i++)
            {
                int bit = (int)bitIndex;
                int pos = bit >> 3;
                int cbit = 1 << (bit & 0x7);
                byte sbit = inMemoryFilter == null
                    ? buffer.get(pos)
                    : inMemoryFilter[pos];

                if ((sbit & cbit) == 0)
                    return false;

                bitIndex += step;
                if (bitIndex >= numBits)
                    bitIndex -= numBits;
            }
            return true;
        }

        public long getNumElements()
        {
            return numElements;
        }

        public long getEstimatedDistinctElements()
        {
            long setBits = 0;
            int numBytes = (int) (((long) numBits + 7L) >>> 3);

            if (inMemoryFilter != null)
            {
                for (int i = 0; i < numBytes; i++)
                    setBits += BYTE_COUNTS[inMemoryFilter[i] & 0xFF];
            }
            else
            {
                for (int i = 0; i < numBytes; i++)
                    setBits += BYTE_COUNTS[buffer.get(i) & 0xFF];
            }

            if (setBits == 0)
                return 0;
            if (setBits >= numBits)
                return Long.MAX_VALUE;

            double estimate =
                -(double)numBits / (double)numHashes
                * Math.log(1.0 - (double)setBits / (double)numBits);
            return Math.round(estimate);
        }

        public void close() throws IOException
        {
            if (inMemoryFilter != null && heapBuildFilename != null && !readOnly)
            {
                int fileLength = inMemoryFilter.length;
                ByteBuffer.wrap(
                    inMemoryFilter,
                    fileLength - NUM_ELEMENTS_OFFSET_FROM_END,
                    Long.BYTES
                ).putLong(numElements);
                File target = new File(heapBuildFilename);
                if (target.exists() && !target.delete())
                    throw new IOException("Could not replace existing file: " + heapBuildFilename);
                try (BufferedOutputStream out = new BufferedOutputStream(
                        new FileOutputStream(target), 1 << 20))
                {
                    out.write(inMemoryFilter);
                }
                heapBuildFilename = null;
                inMemoryFilter = null;
                return;
            }
            if (buffer != null)
            {
                if (!readOnly)
                {
                    int fileLength = (int)file.length();
                    buffer.putLong(
                        fileLength - NUM_ELEMENTS_OFFSET_FROM_END,
                        numElements
                    );
                    buffer.force();
                }

                buffer.clear();
                closeDirectBuffer(buffer);
                if (file != null)
                    file.close();
                buffer = null;
                file = null;
            }
            else if (file != null)
            {
                file.close();
                file = null;
            }
        }

        protected void closeDirectBuffer(ByteBuffer bb)
        {
            if (!bb.isDirect())
                return;

            try
            {
                Method cleanerMethod = bb.getClass().getMethod("cleaner");
                cleanerMethod.setAccessible(true);
                Object cleaner = cleanerMethod.invoke(bb);
                if (cleaner != null)
                {
                    Method cleanMethod = cleaner.getClass().getMethod("clean");
                    cleanMethod.invoke(cleaner);
                }
            }
            catch (Exception e)
            {
                // Best-effort unmapping; the JVM can clean the buffer later.
            }
        }
    }
}
