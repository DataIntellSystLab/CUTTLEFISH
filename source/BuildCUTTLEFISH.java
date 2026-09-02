//Classification by U-mer Tuples of Taxa with Likelihood Estimated Filters Incrementally Structured under a Hierarchy (CUTTLEFISH)
//Checkpoint: COSLOG percentile ln-odds strict-SMH build with exact rolling canonical
//m-mer encoding, unrolled step-2/4/6 minimizer scans, and concise production
//clustering progress output, hash-once shard-grouped Bloom insertion batches,
//and canonical genome-ID ordering for reproducible seeded clustering.

import java.io.*;
import java.io.File;
import java.lang.*;
import java.math.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.nio.*;
import java.nio.file.*;
import java.nio.charset.*;
import java.util.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.*;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.LongAdder;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.ToLongFunction;
import com.dynatrace.hash4j.similarity.*;
import org.apache.datasketches.theta.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import it.unimi.dsi.fastutil.longs.LongIterator;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.SynchronousQueue;

public class BuildCUTTLEFISH
{
	static final class BuildMetrics {
		final long wallStartNanos=System.nanoTime();
		final LongAdder fastaDiscoveryNanos=new LongAdder(), fastaReadingNanos=new LongAdder(), minimizerHashingNanos=new LongAdder(), sketchUpdateNanos=new LongAdder(), tempWriteNanos=new LongAdder(), similarityNanos=new LongAdder(), treeBuildNanos=new LongAdder(), topologyBuildNanos=new LongAdder(), bloomPopulationPassNanos=new LongAdder(), clusteringNanos=new LongAdder(), sketchCapacityNanos=new LongAdder(), bloomAllocationNanos=new LongAdder(), tempReadNanos=new LongAdder(), bloomInsertionNanos=new LongAdder(), bloomEstimationNanos=new LongAdder(), bloomCloseNanos=new LongAdder(), serializationNanos=new LongAdder();
		final LongAdder inputBases=new LongAdder(), uniqueMinimizers=new LongAdder(), similarityCalculations=new LongAdder(), bloomInsertions=new LongAdder(), tempBytesWritten=new LongAdder(), tempBytesRead=new LongAdder();
	}
	static final BuildMetrics METRICS=new BuildMetrics();
	private static final ZeroSimilarityCache ZERO_SIMILARITY_CACHE =
		new ZeroSimilarityCache(256);
	private static final LongAdder ZERO_CACHE_HITS = new LongAdder();
	private static final LongAdder ZERO_CACHE_ADDITIONS = new LongAdder();
	private static final LongAdder ZERO_CACHE_DUPLICATE_ADDITIONS = new LongAdder();
	private static final LongAdder ZERO_CACHE_GROWTH_SKIPS = new LongAdder();
	private static final class ZeroSimilarityCache
	{
		private final LongOpenHashSet[] sets;
		private final Object[] locks;
		private final int mask;
		private final AtomicLong addChecks = new AtomicLong();
		private volatile boolean acceptingNewEntries = true;
		ZeroSimilarityCache(int stripeCount)
		{
			if (Integer.bitCount(stripeCount)!=1)
				throw new IllegalArgumentException("stripeCount must be a power of two");
			sets=new LongOpenHashSet[stripeCount];
			locks=new Object[stripeCount];
			mask=stripeCount-1;
			for (int i=0;i<stripeCount;i++)
			{
				sets[i]=new LongOpenHashSet();
				locks[i]=new Object();
			}
		}
		private int stripe(long key)
		{
			long z=key;
			z^=z>>>33;
			z*=0xff51afd7ed558ccdl;
			z^=z>>>33;
			return ((int)z)&mask;
		}
		boolean contains(long key)
		{
			int stripe=stripe(key);
			synchronized(locks[stripe])
			{
				return sets[stripe].contains(key);
			}
		}
		boolean add(long key)
		{
			long checks=addChecks.incrementAndGet();
			if ((checks&0x3fffL)==0L)
			{
				Runtime runtime=Runtime.getRuntime();
				double used=(double)(runtime.totalMemory()-runtime.freeMemory());
				acceptingNewEntries=used/runtime.maxMemory()<0.85d;
			}
			if (!acceptingNewEntries)
			{
				ZERO_CACHE_GROWTH_SKIPS.increment();
				return false;
			}
			int stripe=stripe(key);
			synchronized(locks[stripe])
			{
				return sets[stripe].add(key);
			}
		}
		long size()
		{
			long total=0L;
			for (int i=0;i<sets.length;i++)
			{
				synchronized(locks[i]) { total+=sets[i].size(); }
			}
			return total;
		}
		int maxStripeSize()
		{
			int maximum=0;
			for (int i=0;i<sets.length;i++)
			{
				synchronized(locks[i]) { maximum=Math.max(maximum,sets[i].size()); }
			}
			return maximum;
		}
		boolean isAcceptingNewEntries() { return acceptingNewEntries; }
	}

	static final int THETA_K=4096;
	private static final int SMH_COMPONENTS=4096, SMH_BITS=8;
	static final long THETA_SEED=9001L;
	private static final SimilarityHashPolicy SMH_POLICY=SimilarityHashing.superMinHash(SMH_COMPONENTS,SMH_BITS);
	private static final double SMH_COLLISION=Math.scalb(1.0,-SMH_BITS);
	private static final LongAdder THETA_POSITIVES=new LongAdder(), SMH_UNCERTAIN=new LongAdder(), DUAL_ZEROS=new LongAdder(), SECOND_PASS_SMH=new LongAdder();
	private static long CLUSTERING_SEED=9001L;
	private static final int CALIBRATION_TARGET_POSITIVES=32;
	private static final ConcurrentHashMap<Long,Boolean> CALIBRATION_PAIRS=new ConcurrentHashMap<Long,Boolean>();
	private static final ConcurrentHashMap<Long,Byte> CALIBRATION_SIMILARITY_CLASS=new ConcurrentHashMap<Long,Byte>();
	private static final class LongList extends AbstractList<Long>{final long[] a;LongList(long[] a){this.a=a;}public Long get(int i){return a[i];}public int size(){return a.length;}}
	static final class GenomeSketch{final CompactSketch theta;final byte[] smh;GenomeSketch(CompactSketch t,byte[] s){theta=t;smh=s;}long cardinality(){return Math.max(0L,Math.round(theta.getEstimate()));}}
	private static GenomeSketch buildGenomeSketch(LongOpenHashSet set){long t=System.nanoTime();long[] a=set.toLongArray();UpdateSketch u=UpdateSketch.builder().setNominalEntries(THETA_K).setSeed(THETA_SEED).build();for(long v:a)u.update(v);SimilarityHasher h=SMH_POLICY.createHasher();ToLongFunction<Long> id=x->x;byte[] smh=h.compute(ElementHashProvider.ofCollection(new LongList(a),id));METRICS.sketchUpdateNanos.add(System.nanoTime()-t);return new GenomeSketch(u.compact(),smh);}
	private static double thetaCosineLog(GenomeSketch a,GenomeSketch b){Intersection ix=SetOperation.builder().setSeed(THETA_SEED).buildIntersection();ix.intersect(a.theta);ix.intersect(b.theta);double in=ix.getResult().getEstimate();if(in<=0)return 0;double ca=a.theta.getEstimate(),cb=b.theta.getEstimate(),den=Math.sqrt(ca*cb);if(den<=0)return 0;double cosine=Math.max(0,Math.min(1,in/den));if(cosine<=0)return 0;double score=1.0/(1.0-Math.log(cosine));return Double.isFinite(score)&&score>0.0?Math.min(1.0,score):0.0;}
	private static double smhJaccard(GenomeSketch a,GenomeSketch b){double f=SMH_POLICY.getFractionOfEqualComponents(a.smh,b.smh);return Math.max(0,Math.min(1,(f-SMH_COLLISION)/(1-SMH_COLLISION)));}
	static long thetaUnionCardinality(List<GenomeSketch> x){Union u=SetOperation.builder().setNominalEntries(THETA_K).setSeed(THETA_SEED).buildUnion();for(GenomeSketch g:x)u.union(g.theta);return Math.max(0L,Math.round(u.getResult().getEstimate()));}

	private static double seconds(LongAdder v){return v.sum()/1_000_000_000.0;}
	private static void printPhaseSummary(){
		long w=METRICS.tempBytesWritten.sum();
		double ratio=w==0?0.0:(double)METRICS.tempBytesRead.sum()/w;
		long treeTotal=METRICS.treeBuildNanos.sum();
		long treeMeasured=METRICS.clusteringNanos.sum()+METRICS.sketchCapacityNanos.sum()+METRICS.bloomAllocationNanos.sum()+METRICS.tempReadNanos.sum()+METRICS.bloomInsertionNanos.sum()+METRICS.bloomEstimationNanos.sum()+METRICS.bloomCloseNanos.sum();
		long treeOther=Math.max(0L,treeTotal-treeMeasured);
		System.out.println("\nBuild phase summary");
		System.out.printf(Locale.ROOT,"\tFASTA discovery             = %.3f s%n",seconds(METRICS.fastaDiscoveryNanos));
		System.out.printf(Locale.ROOT,"\tFASTA reading               = %.3f s%n",seconds(METRICS.fastaReadingNanos));
		System.out.printf(Locale.ROOT,"\tMinimizer hashing           = %.3f s%n",seconds(METRICS.minimizerHashingNanos));
		System.out.printf(Locale.ROOT,"\tTheta+SMH updates         = %.3f s%n",seconds(METRICS.sketchUpdateNanos));
		System.out.printf(Locale.ROOT,"\tTemporary writes            = %.3f s%n",seconds(METRICS.tempWriteNanos));
		System.out.printf(Locale.ROOT,"\tSimilarity sampling         = %.3f s%n",seconds(METRICS.similarityNanos));
		System.out.printf(Locale.ROOT,"\tTree/Bloom construction     = %.3f s%n",seconds(METRICS.treeBuildNanos));
		System.out.printf(Locale.ROOT,"\t  Topology-only build       = %.3f s%n",seconds(METRICS.topologyBuildNanos));
		System.out.printf(Locale.ROOT,"\t  Filter-local Bloom pass   = %.3f s%n",seconds(METRICS.bloomPopulationPassNanos));
		System.out.printf(Locale.ROOT,"\t  Genome clustering         = %.3f s%n",seconds(METRICS.clusteringNanos));
		System.out.printf(Locale.ROOT,"\t  Sketch/capacity estimates = %.3f s%n",seconds(METRICS.sketchCapacityNanos));
		System.out.printf(Locale.ROOT,"\t  Bloom allocation/mapping  = %.3f s%n",seconds(METRICS.bloomAllocationNanos));
		System.out.printf(Locale.ROOT,"\t  Temporary reads           = %.3f s%n",seconds(METRICS.tempReadNanos));
		System.out.printf(Locale.ROOT,"\t  Bloom insertion           = %.3f s%n",seconds(METRICS.bloomInsertionNanos));
		System.out.printf(Locale.ROOT,"\t  Bloom cardinality scans   = %.3f s%n",seconds(METRICS.bloomEstimationNanos));
		System.out.printf(Locale.ROOT,"\t  Bloom flush/close          = %.3f s%n",seconds(METRICS.bloomCloseNanos));
		System.out.printf(Locale.ROOT,"\t  Other tree/Bloom work     = %.3f s%n",treeOther/1_000_000_000.0);
		System.out.printf(Locale.ROOT,"\tSerialization/metadata      = %.3f s%n",seconds(METRICS.serializationNanos));
		System.out.println();
		System.out.printf(Locale.ROOT,"\tInput bases                 = %,d%n",METRICS.inputBases.sum());
		System.out.printf(Locale.ROOT,"\tUnique minimizers           = %,d%n",METRICS.uniqueMinimizers.sum());
		System.out.printf(Locale.ROOT,"\tSimilarity calculations    = %,d%n",METRICS.similarityCalculations.sum());
		System.out.printf(Locale.ROOT,"\tBloom insertions            = %,d%n",METRICS.bloomInsertions.sum());
		System.out.printf(Locale.ROOT,"\tTemporary bytes written    = %,d%n",w);
		System.out.printf(Locale.ROOT,"\tTemporary bytes read       = %,d%n",METRICS.tempBytesRead.sum());
		System.out.printf(Locale.ROOT,"\tTemporary read/write ratio = %.3f%n",ratio);
		System.out.println();
		System.out.println("Exact striped zero-similarity cache");
		System.out.printf(Locale.ROOT,"\tCache hits                  = %,d%n",ZERO_CACHE_HITS.sum());
		System.out.printf(Locale.ROOT,"\tDistinct zero pairs cached  = %,d%n",ZERO_SIMILARITY_CACHE.size());
		System.out.printf(Locale.ROOT,"\tSuccessful additions        = %,d%n",ZERO_CACHE_ADDITIONS.sum());
		System.out.printf(Locale.ROOT,"\tDuplicate additions         = %,d%n",ZERO_CACHE_DUPLICATE_ADDITIONS.sum());
		System.out.printf(Locale.ROOT,"\tGrowth-skipped additions    = %,d%n",ZERO_CACHE_GROWTH_SKIPS.sum());
		System.out.printf(Locale.ROOT,"\tMaximum stripe size         = %,d%n",ZERO_SIMILARITY_CACHE.maxStripeSize());
		System.out.printf(Locale.ROOT,"\tAccepting new entries       = %s%n",ZERO_SIMILARITY_CACHE.isAcceptingNewEntries());
		System.out.printf(Locale.ROOT,"\tTotal measured wall time   = %.3f s%n",(System.nanoTime()-METRICS.wallStartNanos)/1_000_000_000.0);
		System.out.printf(Locale.ROOT,"\tTheta-positive comparisons  = %,d%n",THETA_POSITIVES.sum());
		System.out.printf(Locale.ROOT,"\tUncertain SMH positives     = %,d%n",SMH_UNCERTAIN.sum());
		System.out.printf(Locale.ROOT,"\tConfirmed dual zeros        = %,d%n",DUAL_ZEROS.sum());
		System.out.printf(Locale.ROOT,"\tSecond-pass SMH assignments = %,d%n",SECOND_PASS_SMH.sum());
	}

	private static final String DB_FORMAT="2";
	private static void writeMetadata(String outdir,boolean complete,int k,int step,double fpp,int maxKmersPerSpecies,int genomeCount,long bloomFilterCount,long physicalShardCount) throws IOException
	{
		Properties p=new Properties();
		p.setProperty("cuttlefishDbFormat",DB_FORMAT);
		p.setProperty("cuttlefishSoftwareVersion","AUDIT_REMEDIATED_V4_1");
		p.setProperty("buildComplete",Boolean.toString(complete));
		p.setProperty("buildTimestampUtc",java.time.Instant.now().toString());
		p.setProperty("minimizerAlgorithm","canonical-2bit-exact-v1");
		p.setProperty("minimizerWindowDefinition","k-step-adjacent-mmers-v1");
		p.setProperty("bloomFormat","legacy-header-prime-bits-v1");
		p.setProperty("bloomHashAlgorithm","murmurhash2-32-seeds-0-1-big-endian-long");
		p.setProperty("bloomProbePolicy","double-hash-prime-modulus-zero-step-to-one-v1");
		p.setProperty("bloomShardSelection","seed0-low-bits-power-of-two-v1");
		p.setProperty("thetaNominalEntries",Integer.toString(THETA_K));
		p.setProperty("thetaSeed",Long.toString(THETA_SEED));
		p.setProperty("superMinHashComponents",Integer.toString(SMH_COMPONENTS));
		p.setProperty("superMinHashBits",Integer.toString(SMH_BITS));
		p.setProperty("distinctCounterPolicy","auto-ultraloglog-or-datasketches-hll4");
		p.setProperty("distinctHashMixer","murmur3-fmix64-seed-9e3779b97f4a7c15");
		p.setProperty("clusteringSeed",Long.toString(CLUSTERING_SEED));
		p.setProperty("k",Integer.toString(k));p.setProperty("step",Integer.toString(step));
		p.setProperty("fpp",Double.toString(fpp));p.setProperty("maxKmersPerSpecies",Integer.toString(maxKmersPerSpecies));
		p.setProperty("genomeCount",Integer.toString(genomeCount));p.setProperty("bloomFilterCount",Long.toString(bloomFilterCount));
		p.setProperty("physicalShardCount",Long.toString(physicalShardCount));
		Path target=Paths.get(outdir,"metadata.txt"),tmp=Paths.get(outdir,"metadata.txt.tmp");
		try(BufferedWriter writer=Files.newBufferedWriter(tmp,StandardCharsets.UTF_8)){p.store(writer,"CUTTLEFISH database metadata");}
		try{Files.move(tmp,target,StandardCopyOption.REPLACE_EXISTING,StandardCopyOption.ATOMIC_MOVE);}
		catch(AtomicMoveNotSupportedException e){Files.move(tmp,target,StandardCopyOption.REPLACE_EXISTING);}
	}

	public static void main(String[] args) throws Exception
	{
		long timeZero = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;
		
		int k = 29;
		String dbfile="";
		int samplefirstm = -1;
		int cores = Runtime.getRuntime().availableProcessors()-1;
		int maxKmersPerSpecies = Integer.MAX_VALUE-64_424_509; //3% tolerance
		//maxKmersPerSpecies = 1_000_000; //test
		double fpp = 0.00001d;
		long clusteringSeed = 9001L;
		
		float allram = (float)(Runtime.getRuntime().maxMemory());
		System.out.println();
		System.out.println("Max RAM allocated is "+allram/1000000+" MB");
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		
		if (args.length<=0) {System.out.println("type java BuildCUTTLEFISH h [or help, -h, -help] for instructions"); System.exit(0);}
		for (int t = 0; t < args.length; t++)
{
    switch (args[t])
    {
        case "-h":
        case "--help":
            System.out.println(
                "Usage:\n" +
                "\tjava BuildCUTTLEFISH -f <FASTA file or folder> [options]\n\n" +
                "Required:\n" +
                "\t-f <path>       FASTA file, gzipped FASTA, or genome folder\n\n" +
                "Options:\n" +
                "\t-k <integer>    k-mer length, default 29, min 11, max 37\n" +
                "\t-m <integer>    process only the first m genomes\n" +
                "\t-t <integer>    number of worker threads\n" +
                "\t-e <decimal>    Bloom-filter false-positive rate, default 0.00001\n" +
                "\t--clustering-seed <long> deterministic clustering seed, default 9001\n" +
                "\t-h, --help      print this help"
            );
            return;

        case "-f":
            dbfile = requireValue(args, ++t, "-f");
            break;

        case "-k":
            k = Integer.parseInt(
                requireValue(args, ++t, "-k")
            );
            break;

        case "-m":
            samplefirstm = Integer.parseInt(
                requireValue(args, ++t, "-m")
            );

            if (samplefirstm <= 0)
            {
                System.out.println(
                    "Cannot sample fewer than one genome; all genomes will be used."
                );
                samplefirstm = -1;
            }
            break;

        case "-t":
            cores = Integer.parseInt(
                requireValue(args, ++t, "-t")
            );

            if (cores <= 0)
            {
                System.out.println(
                    "Cannot use fewer than one thread; using all processors minus one."
                );
                cores =
                    Runtime.getRuntime().availableProcessors() - 1;
            }
            break;

        case "--clustering-seed":
            clusteringSeed = Long.parseLong(requireValue(args, ++t, "--clustering-seed"));
            break;

        case "-e":
            fpp = Double.parseDouble(
                requireValue(args, ++t, "-e")
            );

            if (fpp <= 0.0 || fpp >= 1.0)
            {
                System.out.println(
                    "Invalid false-positive rate; using 0.00001."
                );
                fpp = 0.00001d;
            }
            break;

        default:
            throw new IllegalArgumentException(
                "Unknown option: " + args[t] +
                ". Use -h for help."
            );
    }
}
		CLUSTERING_SEED=clusteringSeed;
		if (cores<=0) cores=1;
		if (dbfile.equals("")) {System.out.println("Please specify a genomes' input file in fasta format or a folder with all genomes (f:genomes_fasta_file_or_folder)"); System.exit(0);}
		if (k%2==0) k=k+1;
		if (k<11) {System.out.println("Minimum value of k is 11"); k=11;} 
		if (k>37) {System.out.println("Maximum value of k is 37"); k=37;}
		int step = (int)Math.floor(Math.sqrt(k)); if (step%2==1) step=step-1;
		System.out.println();
		System.out.println("Length of kmers is "+k+" (minimizers are k-"+step+", i.e. "+(k-step)+")");
		System.out.println();
				
		System.out.println("False positive rate is "+fpp);
		System.out.println("Running program on "+cores+" threads");
		System.out.println("Deterministic clustering seed is "+CLUSTERING_SEED);
		System.out.println();

		String datafilerootname = dbfile.substring(Math.max(0,dbfile.lastIndexOf("/")+1));
	    //datafilerootname = dbfile.substring(Math.max(0,dbfile.lastIndexOf("\\")+1));System.out.println("nome di1  "+datafilerootname);
		int indperiod = datafilerootname.lastIndexOf(".");
		if (datafilerootname.indexOf(".gz")>-1 && datafilerootname.indexOf(".gz")<indperiod) indperiod=datafilerootname.indexOf(".gz");
		if (datafilerootname.indexOf(".fna")>-1 && datafilerootname.indexOf(".fna")<indperiod) indperiod=datafilerootname.indexOf(".fna");
		if (datafilerootname.indexOf(".fas")>-1 && datafilerootname.indexOf(".fas")<indperiod) indperiod=datafilerootname.indexOf(".fas");
		if (indperiod!=-1) {datafilerootname = datafilerootname.substring(0,indperiod);}
		//File dbfilefile = new File(dbfile); if (dbfilefile.isDirectory()) {datafilerootname = dbfile.substring(Math.max(0,dbfile.lastIndexOf("/")+1));}
		String outdir = datafilerootname+"_CUTTLEFISHdb/";
		File dbdir = new File(outdir);
		System.gc();
		deleteFileOrFolder(dbdir);
		dbdir.mkdirs();
		writeMetadata(outdir,false,k,step,fpp,maxKmersPerSpecies,0,0,0L);
		
		System.out.println("Processing genomes, extracting/writing their kmer spectrum and calculating Theta 4096 + SuperMinHash 4096x8 sketches");
		startTime = System.currentTimeMillis();
		ArrayList<GenomeSketch> genomesHLL = new ArrayList<GenomeSketch>();
		TreeMap<Integer,String> genomesNames = new TreeMap<Integer,String>();
		processGenomes(dbfile, outdir, samplefirstm, k, step, genomesNames, genomesHLL, maxKmersPerSpecies);
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("\tAvailable genomes ("+genomesHLL.size()+") processed in "+elapsedTime/1000+" s");
		System.out.println();
		
		int minSample = (int)( 1000 + genomesHLL.size()*Math.log(genomesHLL.size())*Math.log(genomesHLL.size()) );
		System.out.println("Calculating "+minSample+" samples (1K+NlogNlogN) of cosine-log similarity matrix via Theta 4096");
		startTime = System.currentTimeMillis();
		ConcurrentHashMap<Long,Double> similarityMatrixMap = new ConcurrentHashMap<Long,Double>();
		long similarityStartNanos=System.nanoTime();
		calculateSampleSimilarityMatrix(genomesHLL, cores, similarityMatrixMap, k-step, minSample);
		METRICS.similarityNanos.add(System.nanoTime()-similarityStartNanos);
		System.gc();
		ArrayList<Integer> genomesIds = new ArrayList<Integer>(genomesNames.keySet());
		// ConcurrentHashMap iteration order is unspecified. Canonical ordering makes
		// every downstream seeded shuffle start from the same genome-ID sequence.
		Collections.sort(genomesIds);
		double [] global_avg_std_similarity = calculateAvgStdSimilarityMatrix(similarityMatrixMap, genomesIds);
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("\tSampled similarity matrix calculated in "+elapsedTime/1000+" s ("+similarityMatrixMap.size()+" non-zero similarities found)");
		System.out.println();
		
		System.out.println("Creating taxonomic clusters' tree topology");
		startTime = System.currentTimeMillis();
		long treeBuildStartNanos=System.nanoTime();
		AtomicInteger clusteredGenomes = new AtomicInteger(0);
		AtomicLong numNodes = new AtomicLong(0);
		FileWriter fw = new FileWriter(outdir+"taxonomicClusterTree.txt");
		BufferedWriter bw = new BufferedWriter(fw);
		BloomPopulationPlan bloomPlan = new BloomPopulationPlan(cores);
		GenomesClustersTree.resetConstructionDepthDiagnostics();
		long topologyStartNanos=System.nanoTime();
		GenomesClustersTree gct = new GenomesClustersTree(genomesIds, similarityMatrixMap, genomesHLL, cores, outdir, clusteredGenomes, bw, numNodes, fpp, maxKmersPerSpecies, k-step, global_avg_std_similarity[0], global_avg_std_similarity[1], bloomPlan);
		long topologyElapsedNanos=System.nanoTime()-topologyStartNanos;
		METRICS.topologyBuildNanos.add(topologyElapsedNanos);
		bw.close(); fw.close();
		System.out.printf(Locale.ROOT,
			"\tTopology complete: %,d/%,d genomes assigned; %,d tree nodes; %,d Bloom filters; %.3f s%n",
			clusteredGenomes.get(),genomesHLL.size(),numNodes.get(),bloomPlan.getFilterCount(),topologyElapsedNanos/1_000_000_000.0);
		System.out.printf(Locale.ROOT,"\tMaximum clustering recursion depth = %,d%n",
			GenomesClustersTree.getMaximumConstructionDepth());
		long populationStartNanos=System.nanoTime();
		bloomPlan.populate(gct, outdir, fpp);
		METRICS.bloomPopulationPassNanos.add(System.nanoTime()-populationStartNanos);
		METRICS.treeBuildNanos.add(System.nanoTime()-treeBuildStartNanos);
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.println("\tTaxonomic tree and Bloom filters built in "+elapsedTime/1000+" s");
		System.out.println("\tTree contains "+numNodes.get()+" nodes");
		System.out.println();
				
		System.out.println("Deleting temporary files, writing taxonomic clusters' tree and ancillary info on to disk..");
		startTime=System.currentTimeMillis();

		long serializationStartNanos=System.nanoTime();
		FileOutputStream fos = new FileOutputStream(new File(outdir+"taxonomicClusterTree.ser"));
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(gct);
		oos.close();
		fos.close();

		fw = new FileWriter(outdir+"info.txt");
		bw = new BufferedWriter(fw);
		bw.write("idx,nkmerminimiz,name");
		bw.newLine();
		for (Map.Entry<Integer,String> entry : genomesNames.entrySet())
		{
			int genomeId = entry.getKey();
			GenomeSketch genomeHLL = genomesHLL.get(genomeId);
			bw.write(genomeId+","+genomeHLL.cardinality()+","+entry.getValue());
			bw.newLine();
			deleteFileOrFolder(new File (outdir+"genome_"+genomeId));
		}
		bw.close();
		fw.close();
		
		writeMetadata(outdir,true,k,step,fpp,maxKmersPerSpecies,genomesIds.size(),bloomPlan.getFilterCount(),bloomPlan.getTotalPhysicalShards());
		
		METRICS.serializationNanos.add(System.nanoTime()-serializationStartNanos);
		elapsedTime = System.currentTimeMillis() - timeZero;
		System.out.println("Total time employed  = "+elapsedTime/1000+" s");
		printPhaseSummary();
	}
		
	public static double calculateSimilarityAndUpdateMatrix(int a,int b,ArrayList<GenomeSketch> genomesHLL,ConcurrentHashMap<Long,Double> similarityMatrixMap,int kappa)throws Exception
	{int i=Math.min(a,b),j=Math.max(a,b);long key=packGenomePair(i,j);Double sim=similarityMatrixMap.get(key);if(sim==null&&ZERO_SIMILARITY_CACHE.contains(key)){ZERO_CACHE_HITS.increment();return 0;}if(sim==null){GenomeSketch x=genomesHLL.get(i),y=genomesHLL.get(j);METRICS.similarityCalculations.increment();double cosine=thetaCosineLog(x,y);if(cosine>0){THETA_POSITIVES.increment();sim=cosine;float all=(float)Runtime.getRuntime().maxMemory(),used=(float)(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory());if(used/all<0.925f)similarityMatrixMap.put(key,sim);}else if(smhJaccard(x,y)>0){SMH_UNCERTAIN.increment();sim=0.0;}else{DUAL_ZEROS.increment();sim=0.0;if(ZERO_SIMILARITY_CACHE.add(key))ZERO_CACHE_ADDITIONS.increment();else if(ZERO_SIMILARITY_CACHE.contains(key))ZERO_CACHE_DUPLICATE_ADDITIONS.increment();}}return sim;}
	private static double calculateSuperMinHashFallback(int a,int b,ArrayList<GenomeSketch> g){return smhJaccard(g.get(a),g.get(b));}
	private static double calculateThetaOnlyAndUpdateMatrix(int a,int b,ArrayList<GenomeSketch> genomesHLL,ConcurrentHashMap<Long,Double> similarityMatrixMap){int i=Math.min(a,b),j=Math.max(a,b);long key=packGenomePair(i,j);Double sim=similarityMatrixMap.get(key);if(sim!=null&&sim>0.0)return sim.doubleValue();if(ZERO_SIMILARITY_CACHE.contains(key)){ZERO_CACHE_HITS.increment();return 0.0;}METRICS.similarityCalculations.increment();double theta=thetaCosineLog(genomesHLL.get(i),genomesHLL.get(j));if(theta>0.0){THETA_POSITIVES.increment();float all=(float)Runtime.getRuntime().maxMemory(),used=(float)(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory());if(used/all<0.925f)similarityMatrixMap.put(key,theta);return theta;}return 0.0;}
	private static double calculateSuperMinHashFallbackAndConfirmZero(int a,int b,ArrayList<GenomeSketch> g){int i=Math.min(a,b),j=Math.max(a,b);long key=packGenomePair(i,j);if(ZERO_SIMILARITY_CACHE.contains(key))return 0.0;double smh=smhJaccard(g.get(i),g.get(j));if(smh>0.0){SMH_UNCERTAIN.increment();return smh;}DUAL_ZEROS.increment();if(ZERO_SIMILARITY_CACHE.add(key))ZERO_CACHE_ADDITIONS.increment();else if(ZERO_SIMILARITY_CACHE.contains(key))ZERO_CACHE_DUPLICATE_ADDITIONS.increment();return 0.0;}

	public static double[] calculateAvgStdSimilarityMatrix(ConcurrentHashMap<Long,Double> similarityMatrixMap, ArrayList<Integer> genomeIds)
	{
		double[] avgstd = new double[2];
		avgstd[0] = -1;
		avgstd[1] = -1;
		if (genomeIds.size()<3) {return avgstd;}
		HashSet<Integer> genomeIdSet = new HashSet<Integer>(genomeIds);
		double avg_similarity = 0;
		ArrayList<Double> std_similarity_vec = new ArrayList<Double>();
		double std_similarity = 0;
		for (Map.Entry<Long,Double> entry : similarityMatrixMap.entrySet())
		{
			long pair = entry.getKey();
			int ind1 = (int)(pair >>> 32);
			int ind2 = (int)pair;
			if ( genomeIdSet.contains(ind1) && genomeIdSet.contains(ind2) )
			{
				double val = entry.getValue();
				avg_similarity += val;
				std_similarity_vec.add(val);
			}
		}
		if (std_similarity_vec.size()<2) {return avgstd;}
		avg_similarity=avg_similarity/std_similarity_vec.size();
		for (int i=0; i<std_similarity_vec.size(); i++) {std_similarity+=(avg_similarity-std_similarity_vec.get(i))*(avg_similarity-std_similarity_vec.get(i));}
		std_similarity = Math.sqrt(std_similarity/(std_similarity_vec.size()-1));
		avgstd[0] = avg_similarity;
		avgstd[1] = std_similarity;
		return avgstd;
	}
	
	public static double pValue(double x, double m, double s)
	{
    		if (!Double.isFinite(x) || !Double.isFinite(m) || !Double.isFinite(s)) {return 1.0;}
		if (s <= 0.0) {return x >= m ? 0.0 : 1.0;}
		NormalDistribution n = new NormalDistribution(m, s);
		return 1.0 - n.cumulativeProbability(x);
	}
	
public static int calculateAdaptiveClusterBound(
		ConcurrentHashMap<Long, Double> similarityMatrixMap,
		ArrayList<Integer> genomeIds,
		double globalAvgSimilarity,
		double globalStdSimilarity)
{
	if (genomeIds == null || genomeIds.isEmpty())
	{
		return 0;
	}

	int n = genomeIds.size();

	if (n == 1)
	{
		return 1;
	}

	double logBound = Math.log1p(n);

	// Safe default if no usable similarity observations are available.
	int maxNumClusters = Math.min(
			n,
			(int) Math.ceil(logBound)
	);

	double[] localAvgStd = calculateAvgStdSimilarityMatrix(
			similarityMatrixMap,
			genomeIds
	);

	if (localAvgStd == null
			|| localAvgStd.length < 2
			|| !Double.isFinite(localAvgStd[0])
			|| !Double.isFinite(localAvgStd[1])
			|| localAvgStd[0] < 0.0
			|| localAvgStd[1] < 0.0
			|| !Double.isFinite(globalAvgSimilarity)
			|| !Double.isFinite(globalStdSimilarity)
			|| globalAvgSimilarity < 0.0
			|| globalStdSimilarity <= 0.0)
	{
		return maxNumClusters;
	}

	double localAvgSimilarity = localAvgStd[0];
	double localStdSimilarity = localAvgStd[1];

	/*
	 * Evidence from a lower local average similarity.
	 */
	double meanDiversityScore = Math.max(
			0.0,
			(globalAvgSimilarity - localAvgSimilarity)
					/ globalStdSimilarity
	);

	/*
	 * Evidence from a wider local similarity distribution.
	 *
	 * localStd == globalStd -> 0
	 * localStd == 2 * globalStd -> 1
	 * localStd == 3 * globalStd -> 2
	 */
	double spreadDiversityScore = Math.max(
			0.0,
			(localStdSimilarity - globalStdSimilarity)
					/ globalStdSimilarity
	);

	double diversityMultiplier =
			1.0
			+ meanDiversityScore
			+ spreadDiversityScore;

	double adaptiveBound =
			diversityMultiplier * logBound;

	if (!Double.isFinite(adaptiveBound))
	{
		return maxNumClusters;
	}

	double sqrtBound = Math.sqrt(n);

	double finalBound = Math.min(
			adaptiveBound,
			sqrtBound
	);

	maxNumClusters = Math.min(
			n,
			(int) Math.ceil(finalBound)
	);

	return Math.min(n, Math.max(2, maxNumClusters));
}
	
	private static final class CalibrationEvaluation {
		final double thetaScore; final byte similarityClass;
		CalibrationEvaluation(double thetaScore,byte similarityClass){this.thetaScore=thetaScore;this.similarityClass=similarityClass;}
	}
	private static CalibrationEvaluation evaluateCalibrationPair(int a,int b,ArrayList<GenomeSketch> sketches,ConcurrentHashMap<Long,Double> map){int i=Math.min(a,b),j=Math.max(a,b);long key=packGenomePair(i,j);Double cached=map.get(key);if(cached!=null&&cached>0.0){CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)2));return new CalibrationEvaluation(cached.doubleValue(),(byte)2);}Byte known=CALIBRATION_SIMILARITY_CLASS.get(key);if(known!=null){return new CalibrationEvaluation(0.0,known.byteValue());}if(ZERO_SIMILARITY_CACHE.contains(key)){CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)1));return new CalibrationEvaluation(0.0,(byte)1);}GenomeSketch x=sketches.get(i),y=sketches.get(j);METRICS.similarityCalculations.increment();double theta=thetaCosineLog(x,y);if(theta>0.0){THETA_POSITIVES.increment();map.put(key,theta);CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)2));return new CalibrationEvaluation(theta,(byte)2);}double smh=smhJaccard(x,y);if(smh>0.0){SMH_UNCERTAIN.increment();CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)2));return new CalibrationEvaluation(0.0,(byte)2);}DUAL_ZEROS.increment();if(ZERO_SIMILARITY_CACHE.add(key))ZERO_CACHE_ADDITIONS.increment();else if(ZERO_SIMILARITY_CACHE.contains(key))ZERO_CACHE_DUPLICATE_ADDITIONS.increment();CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)1));return new CalibrationEvaluation(0.0,(byte)1);}
	private static long mixSeed(long z){z^=z>>>33;z*=0xff51afd7ed558ccdL;z^=z>>>33;z*=0xc4ceb9fe1a85ec53L;return z^(z>>>33);}
	private static long nodeSeed(List<Integer> ids,long tag){long z=CLUSTERING_SEED^tag;for(int id:ids)z=mixSeed(z^id);return z;}
	private static ArrayList<Double> calibrationPositiveValues(Map<Long,Double> map,List<Integer> ids){HashSet<Integer> set=new HashSet<Integer>(ids);ArrayList<Double> v=new ArrayList<Double>();for(long p:CALIBRATION_PAIRS.keySet()){if(set.contains((int)(p>>>32))&&set.contains((int)p)){Double x=map.get(p);if(x!=null&&x>0)v.add(x);}}Collections.sort(v);return v;}
	private static int calibrationPairCount(List<Integer> ids){HashSet<Integer> set=new HashSet<Integer>(ids);int n=0;for(long p:CALIBRATION_PAIRS.keySet())if(set.contains((int)(p>>>32))&&set.contains((int)p))n++;return n;}
	private static long[] calibrationDecisiveCounts(Map<Long,Double> map,List<Integer> ids){HashSet<Integer> set=new HashSet<Integer>(ids);long zero=0,positive=0,uncertain=0,missing=0;for(long p:CALIBRATION_PAIRS.keySet()){if(!set.contains((int)(p>>>32))||!set.contains((int)p))continue;Double theta=map.get(p);if(theta!=null&&theta>0.0)positive++;else if(ZERO_SIMILARITY_CACHE.contains(p))zero++;else{Byte similarityClass=CALIBRATION_SIMILARITY_CLASS.get(p);if(similarityClass!=null&&similarityClass.byteValue()==2)uncertain++;else missing++;}}return new long[]{zero,positive,uncertain,missing};}
	// Parameter-free natural-log odds compression for the first-pass leader ceiling.
	// P is positive Theta only; Z is confirmed dual zero only. SMH-only outcomes
	// are excluded. Exact endpoints are handled explicitly.
	private static double naturalLogOddsGrowthWeight(long zero,long positive){
		if(positive<=0L&&zero<=0L)return 0.0;
		if(positive<=0L)return 0.0;
		if(zero<=0L)return 1.0;
		double logOdds=Math.log((double)positive/(double)zero);
		return 0.5+0.5*(logOdds/(1.0+Math.abs(logOdds)));
	}
	private static int weightedFirstPassLeaderLimit(int n,long zero,long positive){
		if(n<=1)return n;
		double root=Math.sqrt(n);
		double weight=naturalLogOddsGrowthWeight(zero,positive);
		double weighted=root+((double)n-root)*weight;
		return Math.max(1,Math.min(n,(int)Math.ceil(weighted)));
	}
	private static double positivePercentile(ArrayList<Double> sorted,double score){if(score<=0||sorted.isEmpty())return 0;int lo=0,hi=sorted.size();while(lo<hi){int m=(lo+hi)>>>1;if(sorted.get(m)<=score)lo=m+1;else hi=m;}return lo/(double)sorted.size();}
	private static double deterministicUnit(long seed){return ((mixSeed(seed)>>>11)*0x1.0p-53);}
	private static int augmentCalibration(ArrayList<Integer> ids,ArrayList<GenomeSketch> sketches,ConcurrentHashMap<Long,Double> map,int kappa,int budget)throws Exception{if(ids.size()<2||budget<=0)return 0;ArrayList<Long> pairs=new ArrayList<Long>();for(int i=0;i<ids.size();i++)for(int j=i+1;j<ids.size();j++){long p=packGenomePair(ids.get(i),ids.get(j));if(!CALIBRATION_PAIRS.containsKey(p))pairs.add(p);}Collections.shuffle(pairs,new Random(nodeSeed(ids,0x43414c494252L)));int used=0,positive=calibrationPositiveValues(map,ids).size(),target=Math.min(CALIBRATION_TARGET_POSITIVES,ids.size()*(ids.size()-1)/2);for(long p:pairs){if(used>=budget||positive>=target)break;CALIBRATION_PAIRS.put(p,Boolean.TRUE);CalibrationEvaluation e=evaluateCalibrationPair((int)(p>>>32),(int)p,sketches,map);if(e.thetaScore>0.0)positive++;used++;}return used;}

	private static ArrayList<Integer> deterministicRepresentatives(List<Integer> cluster,int count,long seed){ArrayList<Integer>x=new ArrayList<Integer>(cluster);Collections.shuffle(x,new Random(seed));if(x.size()>count)x.subList(count,x.size()).clear();return x;}
	private static ArrayList<ArrayList<Integer>> deterministicSplit(ArrayList<Integer> ids){ArrayList<Integer>x=new ArrayList<Integer>(ids);Collections.shuffle(x,new Random(nodeSeed(ids,0x53504c4954L)));ArrayList<Integer>a=new ArrayList<Integer>(),b=new ArrayList<Integer>();for(int i=0;i<x.size();i++)(i%2==0?a:b).add(x.get(i));ArrayList<ArrayList<Integer>>r=new ArrayList<ArrayList<Integer>>();r.add(a);if(!b.isEmpty())r.add(b);return r;}

	private static ArrayList<ArrayList<Integer>> validatedPartition(List<Integer> source,ArrayList<ArrayList<Integer>> clusters,String stage)
	{HashSet<Integer> expected=new HashSet<Integer>(source),seen=new HashSet<Integer>();int count=0;for(ArrayList<Integer> c:clusters){if(c==null||c.isEmpty())throw new IllegalStateException(stage+": empty cluster");for(int id:c){count++;if(!expected.contains(id))throw new IllegalStateException(stage+": foreign genome "+id);if(!seen.add(id))throw new IllegalStateException(stage+": duplicate genome "+id);}}if(count!=source.size()||seen.size()!=expected.size())throw new IllegalStateException(stage+": partition mismatch source="+source.size()+" assigned="+count+" unique="+seen.size());return clusters;}

	public static ArrayList<ArrayList<Integer>> clusterGenomes(
		ArrayList<Integer> genomesIds,
		ConcurrentHashMap<Long,Double> similarityMatrixMap,
		ArrayList<GenomeSketch> genomesHLL,
		int numThreads,
		int kappa,
		double globalAvgSimilarity,
		double globalStdSimilarity) throws Exception
	{
		final int minClustSize=10;
		if(genomesIds.size()<minClustSize)
		{
			ArrayList<ArrayList<Integer>> result=new ArrayList<ArrayList<Integer>>();
			result.add(new ArrayList<Integer>(genomesIds));
			return validatedPartition(genomesIds,result,"below-minimum");
		}

		ArrayList<Integer> ids=new ArrayList<Integer>(genomesIds);
		Collections.sort(ids);
		long possible=(long)ids.size()*(ids.size()-1L)/2L;
		int budget=(int)Math.min(possible,(long)Math.ceil(Math.pow(Math.log(Math.max(2,ids.size())),3)));
		ArrayList<Double> calibration=calibrationPositiveValues(similarityMatrixMap,ids);
		int added=0;
		if(calibration.size()<Math.min((long)CALIBRATION_TARGET_POSITIVES,possible))
		{
			added=augmentCalibration(ids,genomesHLL,similarityMatrixMap,kappa,budget);
			calibration=calibrationPositiveValues(similarityMatrixMap,ids);
		}
		long[] calibrationCounts=calibrationDecisiveCounts(similarityMatrixMap,ids);
		long calibrationZero=calibrationCounts[0],calibrationPositive=calibrationCounts[1],calibrationUncertain=calibrationCounts[2],calibrationMissing=calibrationCounts[3];
		long decisive=calibrationZero+calibrationPositive;
		long classifiedTotal=decisive+calibrationUncertain+calibrationMissing;
		long calibrationRegistryCount=calibrationPairCount(ids);
		if(classifiedTotal!=calibrationRegistryCount) throw new IllegalStateException("Calibration classification mismatch: registry="+calibrationRegistryCount+" classified="+classifiedTotal);
		double growthWeight=naturalLogOddsGrowthWeight(calibrationZero,calibrationPositive);
		double logOdds=(calibrationPositive>0L&&calibrationZero>0L)?Math.log((double)calibrationPositive/(double)calibrationZero):(calibrationPositive>0L?Double.POSITIVE_INFINITY:(calibrationZero>0L?Double.NEGATIVE_INFINITY:Double.NaN));
		int firstPassLeaderLimit=weightedFirstPassLeaderLimit(ids.size(),calibrationZero,calibrationPositive);

		if(calibration.isEmpty())
		{
			return validatedPartition(ids,deterministicSplit(ids),"zero-calibration-split");
		}

		// Calibration evidence is frozen before leader comparisons.  P uses only positive Theta scores in the map; Z uses only confirmed dual zeros in the zero cache.  SMH-only positives are uncertain and excluded from the weighted-limit denominator.
		ArrayList<ArrayList<Integer>> provisional=new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> deferred=new ArrayList<ArrayList<Integer>>();
		long nodeRandomSeed=nodeSeed(ids,0x50455243454e54L);
		ArrayList<Integer> leaderOrder=new ArrayList<Integer>(ids);
		Collections.shuffle(leaderOrder,new Random(nodeSeed(ids,0x4c45414445524fL)));
		for(int item:leaderOrder)
		{
			boolean placed=false;
			ArrayList<Integer> clusterOrder=new ArrayList<Integer>();
			for(int i=0;i<provisional.size();i++) clusterOrder.add(i);
			Collections.shuffle(clusterOrder,new Random(mixSeed(nodeRandomSeed^item)));
			for(int clusterIndex:clusterOrder)
			{
				ArrayList<Integer> cluster=provisional.get(clusterIndex);
				int representativeCount=(int)Math.sqrt(cluster.size())+1;
				for(int representative:deterministicRepresentatives(
					cluster,representativeCount,
					mixSeed(nodeRandomSeed^item^((long)clusterIndex<<32))))
				{
					double score=calculateSimilarityAndUpdateMatrix(
						item,representative,genomesHLL,similarityMatrixMap,kappa);
					double percentile=positivePercentile(calibration,score);
					double draw=deterministicUnit(
						nodeRandomSeed^((long)item<<32)^representative^clusterIndex);
					if(score>0.0&&draw<percentile)
					{
						cluster.add(item);
						placed=true;
						break;
					}
				}
				if(placed) break;
			}
			if(!placed)
			{
				ArrayList<Integer> cluster=new ArrayList<Integer>();
				cluster.add(item);
				if(provisional.size()<firstPassLeaderLimit) provisional.add(cluster);
				else deferred.add(cluster);
			}
		}


		// Mean/std are intentionally calculated after the leader pass, so the
		// adaptive bound can use both calibration and newly discovered positives.
		int maxNumClusters=calculateAdaptiveClusterBound(
			similarityMatrixMap,ids,globalAvgSimilarity,globalStdSimilarity);
		maxNumClusters=Math.max(3,maxNumClusters);
		maxNumClusters=Math.min(maxNumClusters,provisional.size());

		// Seeded order plus stable size sort: equal-size ties retain seeded order.
		Collections.shuffle(provisional,new Random(nodeSeed(ids,0x52455441494eL)));
		Collections.sort(provisional,new Comparator<ArrayList<Integer>>() {
			@Override public int compare(ArrayList<Integer> a,ArrayList<Integer> b)
			{return Integer.compare(b.size(),a.size());}
		});

		ArrayList<ArrayList<Integer>> retained=new ArrayList<ArrayList<Integer>>();
		int retainedCount=Math.min(maxNumClusters,provisional.size());
		for(int i=0;i<retainedCount;i++) retained.add(provisional.get(i));

		// First-pass overflow genomes and provisional leaders beyond the final
		// adaptive maximum share the same bounded second-pass assignment.
		for(int i=retainedCount;i<provisional.size();i++) deferred.add(provisional.get(i));
		for(ArrayList<Integer> loose:deferred)
		{
			int item=loose.get(0);
			int bestThetaIndex=-1;
			double bestThetaScore=0.0;
			ArrayList<ArrayList<Integer>> representativePanels=new ArrayList<ArrayList<Integer>>(retained.size());
			for(int clusterIndex=0;clusterIndex<retained.size();clusterIndex++)
			{
				ArrayList<Integer> target=retained.get(clusterIndex);
				int representativeCount=(int)Math.sqrt(target.size())+1;
				ArrayList<Integer> panel=deterministicRepresentatives(
					target,representativeCount,
					mixSeed(nodeRandomSeed^0x5345434f4e44L^item^((long)clusterIndex<<32)));
				representativePanels.add(panel);
				for(int representative:panel)
				{
					double thetaScore=calculateThetaOnlyAndUpdateMatrix(
						item,representative,genomesHLL,similarityMatrixMap);
					if(thetaScore>bestThetaScore)
					{
						bestThetaScore=thetaScore;
						bestThetaIndex=clusterIndex;
					}
				}
			}

			if(bestThetaIndex>=0)
			{
				retained.get(bestThetaIndex).addAll(loose);
				continue;
			}

			// Strict fallback: SuperMinHash is evaluated only after every bounded
			// Theta comparison against every retained cluster has returned zero.
			int bestSmhIndex=-1;
			double bestSmhScore=0.0;
			for(int clusterIndex=0;clusterIndex<representativePanels.size();clusterIndex++)
			{
				for(int representative:representativePanels.get(clusterIndex))
				{
					double smhScore=calculateSuperMinHashFallbackAndConfirmZero(
						item,representative,genomesHLL);
					if(smhScore>bestSmhScore)
					{
						bestSmhScore=smhScore;
						bestSmhIndex=clusterIndex;
					}
				}
			}

			if(bestSmhIndex>=0)
			{
				retained.get(bestSmhIndex).addAll(loose);
				SECOND_PASS_SMH.increment();
			}
			else
			{
				int fallback=(int)Math.floorMod(
					mixSeed(nodeRandomSeed^0x46414c4cL^item),retained.size());
				retained.get(fallback).addAll(loose);
			}
		}

		// Pool every undersized retained cluster into one residual cluster.
		ArrayList<ArrayList<Integer>> consolidated=new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> residual=new ArrayList<Integer>();
		for(ArrayList<Integer> cluster:retained)
		{
			if(cluster.size()<minClustSize) residual.addAll(cluster);
			else consolidated.add(cluster);
		}
		if(!residual.isEmpty()) consolidated.add(residual);


		if(consolidated.size()==1&&consolidated.get(0).size()>=minClustSize)
			return validatedPartition(ids,deterministicSplit(consolidated.get(0)),"single-cluster-split");
		return validatedPartition(ids,consolidated,"consolidated");
	}

	private static final int GENOME_IO_BUFFER_SIZE = 1 << 20;
	// Reused for both heap-serial and mapped-parallel population. This replaces
	// per-genome batch allocation and amortizes temporary-file decoding overhead.
	static final int BLOOM_POPULATION_BATCH_LONGS = 1_000_000;

	private static void longToBytes(long value, byte[] bytes)
	{
		bytes[0] = (byte)(value >>> 56); bytes[1] = (byte)(value >>> 48);
		bytes[2] = (byte)(value >>> 40); bytes[3] = (byte)(value >>> 32);
		bytes[4] = (byte)(value >>> 24); bytes[5] = (byte)(value >>> 16);
		bytes[6] = (byte)(value >>> 8);  bytes[7] = (byte)value;
	}

	private static void retainSuffix(StringBuilder sequence, int length)
	{
		if (sequence.length() > length) sequence.delete(0, sequence.length() - length);
	}
	private static void flushSequenceChunk(StringBuilder sequence,int bigK,int step,long[][] pp,int[] nc,int[] ncr,LongOpenHashSet minimizers)
	{
		if(sequence.length()<bigK)return;
		hashSequenceInto(sequence,bigK,step,pp,nc,ncr,minimizers);
		retainSuffix(sequence,bigK-1);
	}
	private static boolean isGenomeFasta(Path path)
	{
		if (!Files.isRegularFile(path)) return false;
		String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
		if (name.contains("_cds_from_genomic")
			|| name.contains("_rna_from_genomic")
			|| name.contains("_protein")
			|| name.contains("translated_cds")
			|| name.endsWith(".faa")
			|| name.endsWith(".faa.gz"))
		{
			return false;
		}
		return name.endsWith(".fna")
			|| name.endsWith(".fna.gz")
			|| name.endsWith(".fa")
			|| name.endsWith(".fa.gz")
			|| name.endsWith(".fasta")
			|| name.endsWith(".fasta.gz")
			|| name.endsWith(".fas")
			|| name.endsWith(".fas.gz");
	}
	private static ArrayList<Path> findGenomeFastas(Path root) throws IOException
	{
		ArrayList<Path> genomeFiles = new ArrayList<Path>();
		try (java.util.stream.Stream<Path> paths = Files.walk(root))
		{
			paths.filter(BuildCUTTLEFISH::isGenomeFasta).forEach(genomeFiles::add);
		}
		genomeFiles.sort(Comparator.comparing(
			path -> root.relativize(path).toString().replace(File.separatorChar, '/')));
		return genomeFiles;
	}
	private static void appendSequenceLine(StringBuilder sequence, String line)
	{
		long readStartNanos=System.nanoTime(); int before=sequence.length();
		int length = line.length();
		int firstWhitespace = -1;
		for (int i = 0; i < length; i++)
		{
			if (Character.isWhitespace(line.charAt(i)))
			{
				firstWhitespace = i;
				break;
			}
		}
		if (firstWhitespace < 0)
		{
			sequence.append(line); METRICS.inputBases.add(sequence.length()-before); METRICS.fastaReadingNanos.add(System.nanoTime()-readStartNanos);
			return;
		}
		sequence.append(line, 0, firstWhitespace);
		for (int i = firstWhitespace + 1; i < length; i++)
		{
			char nucleotide = line.charAt(i);
			if (!Character.isWhitespace(nucleotide)) sequence.append(nucleotide);
		}
	}

	// Binary temporary format: long count, followed by count long minimizers.
	public static void writeLongSet(LongOpenHashSet values, String pathFile) throws IOException
	{
		long metricStartNanos=System.nanoTime();
		try (DataOutputStream out = new DataOutputStream(new BufferedOutputStream(
			new FileOutputStream(pathFile), GENOME_IO_BUFFER_SIZE)))
		{
			out.writeLong(values.size());
			LongIterator iterator = values.iterator();
			while (iterator.hasNext()) out.writeLong(iterator.nextLong());
		}
		METRICS.tempWriteNanos.add(System.nanoTime()-metricStartNanos);
		METRICS.tempBytesWritten.add(8L+8L*values.size());
		METRICS.uniqueMinimizers.add(values.size());
	}

	public static long readLongSetSize(String filePath) throws IOException
	{
		long metricStartNanos=System.nanoTime();
		try (DataInputStream in = new DataInputStream(new BufferedInputStream(
			new FileInputStream(filePath), GENOME_IO_BUFFER_SIZE)))
		{
			long count = in.readLong();
			if (count < 0) throw new IOException("Invalid minimizer count in " + filePath);
			METRICS.tempBytesRead.add(8); METRICS.tempReadNanos.add(System.nanoTime()-metricStartNanos);
			return count;
		}
	}

	public static long addLongFileToBloom(String filePath, ShardedBloomFilter bf) throws IOException
	{
		try (DataInputStream in = new DataInputStream(new BufferedInputStream(
			new FileInputStream(filePath), GENOME_IO_BUFFER_SIZE)))
		{
			long readStartNanos=System.nanoTime();
			long count = in.readLong();
			METRICS.tempReadNanos.add(System.nanoTime()-readStartNanos);
			if (count < 0) throw new IOException("Invalid minimizer count in " + filePath);
			METRICS.tempBytesRead.add(8L+8L*count);
			long[] batch = new long[(int)Math.min((long)BLOOM_POPULATION_BATCH_LONGS,Math.max(1L,count))];
			long remaining=count;
			while (remaining>0)
			{
				int batchSize=(int)Math.min((long)batch.length,remaining);
				readStartNanos=System.nanoTime();
				for (int i=0;i<batchSize;i++) batch[i]=in.readLong();
				METRICS.tempReadNanos.add(System.nanoTime()-readStartNanos);
				long insertionStartNanos=System.nanoTime();
				bf.addBatchGroupedByShard(batch,batchSize);
				METRICS.bloomInsertionNanos.add(System.nanoTime()-insertionStartNanos);
				METRICS.bloomInsertions.add(batchSize);
				remaining-=batchSize;
			}
			return count;
		}
	}

	public static long readLongSetSizeForPlanning(String filePath) throws IOException
	{
		try (RandomAccessFile in = new RandomAccessFile(filePath, "r"))
		{
			long count = in.readLong();
			if (count < 0) throw new IOException("Invalid minimizer count in " + filePath);
			return count;
		}
	}

	static final class BloomPopulationProgress
	{
		private final long expectedInsertions;
		private final long startNanos = System.nanoTime();
		private long completedInsertions = 0L;
		private int nextPercent = 5;

		BloomPopulationProgress(long expectedInsertions)
		{
			this.expectedInsertions = Math.max(0L, expectedInsertions);
		}

		void add(long insertions)
		{
			completedInsertions += insertions;
			while (nextPercent <= 100 &&
				(expectedInsertions == 0L || completedInsertions * 100L >= expectedInsertions * nextPercent))
			{
				printProgress(nextPercent);
				nextPercent += 5;
			}
		}

		void finish()
		{
			if (nextPercent <= 100) printProgress(100);
			if (expectedInsertions != completedInsertions)
				System.out.println("\tWarning: expected "+expectedInsertions+
					" Bloom insertions but completed "+completedInsertions);
		}

		private void printProgress(int percent)
		{
			double elapsedSeconds = (System.nanoTime()-startNanos)/1_000_000_000.0;
			double throughput = elapsedSeconds <= 0.0 ? 0.0 : completedInsertions/elapsedSeconds;
			double remainingSeconds = throughput <= 0.0 ? Double.NaN :
				Math.max(0L, expectedInsertions-completedInsertions)/throughput;
			String eta = Double.isFinite(remainingSeconds) ?
				formatDuration(remainingSeconds) : "unknown";
			System.out.printf(Locale.ROOT,
				"\t%3d%% Bloom population completed: %,d/%,d insertions; elapsed %s; %.3f M/s; ETA %s%n",
				percent, Math.min(completedInsertions,expectedInsertions), expectedInsertions,
				formatDuration(elapsedSeconds), throughput/1_000_000.0, eta);
		}

		private static String formatDuration(double seconds)
		{
			long total = Math.max(0L, Math.round(seconds));
			long hours = total/3600L;
			long minutes = (total%3600L)/60L;
			long secs = total%60L;
			return hours > 0 ? String.format(Locale.ROOT,"%d:%02d:%02d",hours,minutes,secs) :
				String.format(Locale.ROOT,"%d:%02d",minutes,secs);
		}
	}

	public static long addLongFileToBloom(String filePath, ShardedBloomFilter bf,
		BloomPopulationProgress progress, ExecutorService bloomExecutor,
		long[] batch, byte[] byteBuffer, LongAdder backendInsertionNanos,
		LongAdder backendBatches) throws Exception
	{
		try (DataInputStream in = new DataInputStream(new BufferedInputStream(
			new FileInputStream(filePath), GENOME_IO_BUFFER_SIZE)))
		{
			long readStartNanos=System.nanoTime();
			long count=in.readLong();
			METRICS.tempReadNanos.add(System.nanoTime()-readStartNanos);
			if (count < 0) throw new IOException("Invalid minimizer count in " + filePath);
			METRICS.tempBytesRead.add(8L+8L*count);
			long remaining=count;
			while (remaining>0)
			{
				int batchSize=(int)Math.min((long)batch.length,remaining);
				int bytesNeeded=Math.multiplyExact(batchSize,Long.BYTES);
				readStartNanos=System.nanoTime();
				in.readFully(byteBuffer,0,bytesNeeded);
				for (int i=0,offset=0;i<batchSize;i++,offset+=Long.BYTES)
				{
					batch[i]=((long)(byteBuffer[offset]&0xff)<<56)
						|((long)(byteBuffer[offset+1]&0xff)<<48)
						|((long)(byteBuffer[offset+2]&0xff)<<40)
						|((long)(byteBuffer[offset+3]&0xff)<<32)
						|((long)(byteBuffer[offset+4]&0xff)<<24)
						|((long)(byteBuffer[offset+5]&0xff)<<16)
						|((long)(byteBuffer[offset+6]&0xff)<<8)
						|((long)(byteBuffer[offset+7]&0xff));
				}
				METRICS.tempReadNanos.add(System.nanoTime()-readStartNanos);
				long insertionStartNanos=System.nanoTime();
				if (bf.getNumShards()>1)
					bf.addBatchMappedParallel(batch,batchSize,bloomExecutor);
				else
					bf.addBatchGroupedByShard(batch,batchSize);
				long insertionNanos=System.nanoTime()-insertionStartNanos;
				METRICS.bloomInsertionNanos.add(insertionNanos);
				backendInsertionNanos.add(insertionNanos);
				backendBatches.increment();
				METRICS.bloomInsertions.add(batchSize);
				if (progress != null) progress.add(batchSize);
				remaining-=batchSize;
			}
			return count;
		}
	}

	public static void processGenomes(String dbfile, String outdir, int samplefirstm, int k, int step, TreeMap<Integer,String> genomesNames, ArrayList<GenomeSketch> genomesHLL, int maxKmersPerSpecies) throws Exception
	{
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		long [][] pp = precomputedPowsForNucs(4,31);
		int[] ctbf = new int[65_536];
		for (int i=0; i<ctbf.length; i++) {ctbf[i]=-1;}
			ctbf['a']=0; ctbf['A']=0;
			ctbf['c']=1; ctbf['C']=1;
			ctbf['g']=2; ctbf['G']=2;
			ctbf['t']=3; ctbf['T']=3;
			ctbf['u']=3; ctbf['U']=3;
		final int[] nc = ctbf;
		int[] ctbr = new int[65_536];
		for (int i=0; i<ctbr.length; i++) {ctbr[i]=-1;}
			ctbr['a']=3; ctbr['A']=3;
			ctbr['c']=2; ctbr['C']=2;
			ctbr['g']=1; ctbr['G']=1;
			ctbr['t']=0; ctbr['T']=0;
			ctbr['u']=0; ctbr['U']=0;
		final int[] ncr = ctbr;
		File dbfilefile = new File(dbfile);
		if (!dbfilefile.isDirectory())
		{
			BufferedReader r = new BufferedReader(new FileReader(dbfile),DEFAULT_BUFFER_SIZE);
			if(dbfile.endsWith(".gz")) {r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(dbfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
			String header = "";
			int i=0;
			while (true)
			{
				if (samplefirstm>0 && i>=samplefirstm) {break;}
				if (header==null) {break;}
				while (!header.startsWith(">")) {header = r.readLine(); if (header==null) {break;}}
				StringBuilder sequence = new StringBuilder(1 << 20); 
				//String lastkmer = "";
				LongOpenHashSet mvm = new LongOpenHashSet();
				
				String sq = r.readLine();
				while (sq!=null && !sq.startsWith(">"))
				{
					appendSequenceLine(sequence, sq);
					if (sequence.length()>Math.min(10000000,maxKmersPerSpecies))
					{
						flushSequenceChunk(sequence,k,step,pp,nc,ncr,mvm);
					}
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					if (mvm.size()>maxKmersPerSpecies || usedram/allram>0.92f)
					{
						System.out.println("\tWarning: maximum RAM ("+allram/1000000+" MB) or maximum number of k-mers for a single genome ("+maxKmersPerSpecies+") reached; the genome may be split into two or more");
						genomesHLL.add(buildGenomeSketch(mvm)); writeLongSet(mvm,outdir+"genome_"+i); genomesNames.put(i,header.replaceAll(","," - ")); i++; 
						mvm = new LongOpenHashSet();
						System.gc();
						//break;
					}
					sq = r.readLine();
				}
				if (sequence.length()>=k)
				{
					hashSequenceInto(sequence,k,step,pp,nc,ncr,mvm);
				}
				if (!mvm.isEmpty()) {genomesHLL.add(buildGenomeSketch(mvm)); writeLongSet(mvm,outdir+"genome_"+i); genomesNames.put(i,header.replaceAll(","," - ")); i++;}
				header = sq;
				if (i%1000==0)
				{
					//System.gc();
					elapsedTime = System.currentTimeMillis() - startTime;
					System.out.print("\t"+genomesHLL.size()+" genomes processed in "+elapsedTime/1000+" seconds");
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					System.out.println(" - "+Math.round((double)usedram/1000000)+" MB RAM used ("+Math.round(100*(double)usedram/allram)+"%)");
				}
			}
			r.close();
		}
		if (dbfilefile.isDirectory())
		{
			long discoveryStartNanos=System.nanoTime();
			ArrayList<Path> genomeFiles = findGenomeFastas(dbfilefile.toPath());
			METRICS.fastaDiscoveryNanos.add(System.nanoTime()-discoveryStartNanos);
			if (genomeFiles.isEmpty())
				throw new IOException("No supported genomic FASTA files found under: " + dbfile);
			System.out.println("\tFound " + genomeFiles.size() + " genomic FASTA files recursively");
			int i=0;
			for (int c=0; c<genomeFiles.size(); c++)
			{
				if (samplefirstm>0 && i>=samplefirstm) {break;}
				Path genomePath = genomeFiles.get(c);
				String species = genomePath.toString();
				String displayPath = dbfilefile.toPath().relativize(genomePath).toString();
				BufferedReader rd = new BufferedReader(new FileReader(species),DEFAULT_BUFFER_SIZE);
				if(species.toLowerCase(Locale.ROOT).endsWith(".gz")) {rd=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(species),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);}
				String l = rd.readLine();
				while(l!=null && !l.startsWith(">")) {l=rd.readLine();}
				if (l==null) {System.out.println("Empty or missing FASTA header for file "+displayPath+"; skipping it.");}
				if (l!=null)
				{
					String header = l;
					StringBuilder sequence = new StringBuilder(1 << 20); 
					//String lastkmer = "";
					LongOpenHashSet mvm = new LongOpenHashSet();
					
					String sq = rd.readLine();
					while (sq!=null)
					{
						if (!sq.startsWith(">"))
						{
							appendSequenceLine(sequence, sq);
							if (sequence.length()>Math.min(10000000,maxKmersPerSpecies))
							{
								flushSequenceChunk(sequence,k,step,pp,nc,ncr,mvm);
							}
							usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
							if (mvm.size()>maxKmersPerSpecies || usedram/allram>0.92f)
							{
								System.out.println("\tWarning: maximum RAM ("+allram/1000000+" MB) or maximum number of k-mers for a single genome ("+maxKmersPerSpecies+") reached; the genome may be split into two or more"); 
								//System.out.println(i+" "+header.substring(0,25));
								genomesHLL.add(buildGenomeSketch(mvm)); writeLongSet(mvm,outdir+"genome_"+i); genomesNames.put(i,header.replaceAll(","," - ")); i++; 
								mvm = new LongOpenHashSet();
								System.gc(); 
								//break;
							}
						}
						sq = rd.readLine();
					}
					if (sequence.length()>=k)
					{
						hashSequenceInto(sequence,k,step,pp,nc,ncr,mvm);
					}
					if (!mvm.isEmpty()) {genomesHLL.add(buildGenomeSketch(mvm)); writeLongSet(mvm,outdir+"genome_"+i); genomesNames.put(i,header.replaceAll(","," - ")); i++;}
				}
				rd.close();
				if (c%1000==0)
				{
					//System.gc();
					elapsedTime = System.currentTimeMillis() - startTime;
					System.out.print("\t"+genomesHLL.size()+" genomes processed in "+elapsedTime/1000+" seconds");
					usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
					System.out.println(" - "+Math.round((double)usedram/1000000)+" MB RAM used ("+Math.round(100*(double)usedram/allram)+"%)");
				}
			}
		}
	}
	
	private static long packGenomePair(int a, int b)
	{
		int low = Math.min(a, b);
		int high = Math.max(a, b);
		return ((long)low << 32) | Integer.toUnsignedLong(high);
	}

	private static LongOpenHashSet generateUniqueGenomePairs(int numGenomes,int minSample,int numThreads)
	{long max=(long)numGenomes*(numGenomes-1L)/2L;int target=(int)Math.min(max,(long)minSample);LongOpenHashSet out=new LongOpenHashSet(target);if(max<=target){for(int a=0;a<numGenomes;a++)for(int b=a+1;b<numGenomes;b++)out.add(packGenomePair(a,b));return out;}Random r=new Random(mixSeed(CLUSTERING_SEED^0x494e495449414cL));while(out.size()<target){int a=r.nextInt(numGenomes),b=r.nextInt(numGenomes);if(a!=b)out.add(packGenomePair(a,b));}return out;
	}

	public static void calculateSampleSimilarityMatrix(ArrayList<GenomeSketch> genomesHLL, int numThreads, ConcurrentHashMap<Long,Double> similarityMatrixMap, int kappa, int minSample) throws Exception
	{
		if (genomesHLL.size()<10) {System.out.println("Less than ten genomes?!? Not worth it. Exiting program"); System.exit(0);}
		float allram = (float)(Runtime.getRuntime().maxMemory());
		long startTime = System.currentTimeMillis();
		long maxUniquePairs = (long)genomesHLL.size() * (genomesHLL.size()-1L) / 2L;
		long requestedWithExtra = (long)minSample;
		LongOpenHashSet uniquePairs = generateUniqueGenomePairs(genomesHLL.size(),minSample,numThreads);
		long[] pairs = uniquePairs.toLongArray();
		int numTasks = pairs.length;
		int requiredCompletions = numTasks;
		boolean exhaustive = maxUniquePairs <= minSample;
		System.out.println("	Unique pairs requested = "+Math.min(requestedWithExtra,maxUniquePairs)+
			"; submitted = "+numTasks+"; exhaustive = "+exhaustive);
		AtomicInteger tasksCompleted = new AtomicInteger(0);
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		LinkedList<Future> futList = new LinkedList<Future>();
		int lastProgressDecile = -1;
		for (int index=0; index<numTasks; index++)
		{
			final long packedPair = pairs[index];
			final int ind1 = (int)(packedPair >>> 32);
			final int ind2 = (int)packedPair;
			final GenomeSketch gen1 = genomesHLL.get(ind1);
			final GenomeSketch gen2 = genomesHLL.get(ind2);
			Future fut = executor.submit(() ->
			{
				long key = packedPair;
				CALIBRATION_PAIRS.put(key,Boolean.TRUE);
				METRICS.similarityCalculations.increment();
				double cosine=thetaCosineLog(gen1,gen2);
				if(cosine>0){CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)2));similarityMatrixMap.put(key,cosine);THETA_POSITIVES.increment();}
				else{double smh=smhJaccard(gen1,gen2);if(smh>0){CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)2));SMH_UNCERTAIN.increment();}else{CALIBRATION_SIMILARITY_CLASS.put(key,Byte.valueOf((byte)1));DUAL_ZEROS.increment();if(ZERO_SIMILARITY_CACHE.add(key))ZERO_CACHE_ADDITIONS.increment();else if(ZERO_SIMILARITY_CACHE.contains(key))ZERO_CACHE_DUPLICATE_ADDITIONS.increment();}}
				tasksCompleted.incrementAndGet();
			});
			futList.add(fut);
			float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
			if (usedram / allram > 0.9d || futList.size() > numThreads)
			{
				Iterator<Future> iterator = futList.iterator();
				while (iterator.hasNext())
				{
					Future completed = iterator.next();
					if (completed.isDone())
					{
						iterator.remove();
						try {completed.get();} catch (Exception e) {e.printStackTrace();}
					}
				}
			}
			int completed = tasksCompleted.get();
			int progressPercent = Math.min(100,
				(int)(100L*completed/Math.max(1,requiredCompletions)));
			int progressDecile = progressPercent/10;
			if (progressDecile > lastProgressDecile)
			{
				lastProgressDecile = progressDecile;
				long elapsedTime = System.currentTimeMillis()-startTime;
				System.out.print("\t"+progressPercent+"% similarity target completed ("+
					completed+"/"+requiredCompletions+"; "+numTasks+" unique tasks submitted)");
				System.out.print(" - "+elapsedTime/1000+"s");
				System.out.println(" - "+Math.round((double)usedram/1000000)+" MB RAM used ("+
					Math.round(100*(double)usedram/allram)+"%)");
			}
		}
		long deadline = System.nanoTime() + TimeUnit.MINUTES.toNanos(42);
		while (tasksCompleted.get() < requiredCompletions && System.nanoTime() < deadline)
		{
			Thread.sleep(10);
			int completed = tasksCompleted.get();
			int progressPercent = Math.min(100,
				(int)(100L*completed/Math.max(1,requiredCompletions)));
			int progressDecile = progressPercent/10;
			if (progressDecile > lastProgressDecile)
			{
				lastProgressDecile = progressDecile;
				long elapsedTime = System.currentTimeMillis()-startTime;
				float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
				System.out.print("\t"+progressPercent+"% similarity target completed ("+
					completed+"/"+requiredCompletions+"; "+numTasks+" unique tasks submitted)");
				System.out.print(" - "+elapsedTime/1000+"s");
				System.out.println(" - "+Math.round((double)usedram/1000000)+" MB RAM used ("+
					Math.round(100*(double)usedram/allram)+"%)");
			}
		}
		if (tasksCompleted.get() >= requiredCompletions && lastProgressDecile < 10)
		{
			long elapsedTime = System.currentTimeMillis()-startTime;
			System.out.println("\t100% similarity target completed ("+tasksCompleted.get()+"/"+
				requiredCompletions+"; "+numTasks+" unique tasks submitted) - "+elapsedTime/1000+"s");
		}
		executor.shutdownNow();
		if (tasksCompleted.get() < requiredCompletions)
		{
			System.out.println("		Warning: the minimal unique sample of similarities could not be reached");
		}
	}

	public static long [][] precomputedPowsForNucs(int p, int k)
	{
		long [][] ppn = new long [p][];
		for (int w=0; w<p; w++)
		{
			long [] pp = new long [k];
			for (int i=0; i<k; i++)
			{
				long pow = 1;
				for (int j=0; j<i; j++)
					pow = pow*p;
				pp[i] = pow;
			}
			for (int i=0; i<k; i++) pp[i]=pp[i]*w;
			ppn[w]=pp;
		}
		return ppn;
	}
	public static long polynomialHash(int[] str, int start, int stop, long [][] pp)
	{
		long hash_val = 0;
		int k = stop - start;
		for (int i=0; i<k; i++)
		{
			hash_val = hash_val + pp[str[(stop-1)-i]][i];
		}
		return hash_val;
	}
	public static long nextHash(long [][] pp, long prev_hash, int k, int first, int next)
	{
		long nextH = ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + ( prev_hash - pp[first][k-1] ) + next;
		return nextH;
	}
	// Exact collision-free rolling canonical m-mer encoding followed by a
	// fixed-width direct minimum scan. Under CUTTLEFISH's odd k in [11,37]
	// and even step=floor(sqrt(k)), step is 2, 4, or 6 and m=bigK-step <=31.
	private static long minimumValid2(long a, long b)
	{
		if (a==Long.MIN_VALUE) return b;
		if (b==Long.MIN_VALUE) return a;
		return a<b ? a : b;
	}
	private static long minimumValid4(long a, long b, long c, long d)
	{
		return minimumValid2(minimumValid2(a,b),minimumValid2(c,d));
	}
	private static long minimumValid6(long a, long b, long c, long d, long e, long f)
	{
		return minimumValid2(minimumValid4(a,b,c,d),minimumValid2(e,f));
	}
	private static long[] exactCanonicalMers(CharSequence sequence, int m, int[] nc, int[] ncr)
	{
		int numHashes=sequence.length()-m+1;
		long[] canonical=new long[numHashes];
		Arrays.fill(canonical,Long.MIN_VALUE);
		long mask=(1L<<(2*m))-1L;
		int highShift=2*(m-1);
		long forward=0L, reverse=0L;
		int validRun=0;
		for (int position=0; position<sequence.length(); position++)
		{
			char nucleotide=sequence.charAt(position);
			int code=nucleotide<nc.length ? nc[nucleotide] : -1;
			if (code<0)
			{
				forward=0L;
				reverse=0L;
				validRun=0;
			}
			else
			{
				forward=((forward<<2)|code)&mask;
				reverse=(reverse>>>2)|((long)ncr[nucleotide]<<highShift);
				validRun++;
			}
			if (position>=m-1 && validRun>=m)
			{
				int index=position-m+1;
				canonical[index]=forward<reverse ? forward : reverse;
			}
		}
		return canonical;
	}
	private static void addExactUnrolledMinimizers(long[] canonical, int step, LongOpenHashSet minimizers)
	{
		int windows=canonical.length-step+1;
		if (step==2)
		{
			for (int i=0;i<windows;i++)
			{
				long minimum=minimumValid2(canonical[i],canonical[i+1]);
				if (minimum!=Long.MIN_VALUE) minimizers.add(minimum);
			}
		}
		else if (step==4)
		{
			for (int i=0;i<windows;i++)
			{
				long minimum=minimumValid4(canonical[i],canonical[i+1],canonical[i+2],canonical[i+3]);
				if (minimum!=Long.MIN_VALUE) minimizers.add(minimum);
			}
		}
		else if (step==6)
		{
			for (int i=0;i<windows;i++)
			{
				long minimum=minimumValid6(canonical[i],canonical[i+1],canonical[i+2],canonical[i+3],canonical[i+4],canonical[i+5]);
				if (minimum!=Long.MIN_VALUE) minimizers.add(minimum);
			}
		}
		else
		{
			// Defensive exact fallback; current CUTTLEFISH parameterization never uses it.
			for (int i=0;i<windows;i++)
			{
				long minimum=Long.MIN_VALUE;
				for (int j=0;j<step;j++) minimum=minimumValid2(minimum,canonical[i+j]);
				if (minimum!=Long.MIN_VALUE) minimizers.add(minimum);
			}
		}
	}
	public static void hashSequenceInto(
		CharSequence fwd,
		int bigK,
		int step,
		long [][] pp,
		int[] nc,
		int[] ncr,
		LongOpenHashSet minimizers)
	{
		long hashingStartNanos=System.nanoTime();
		if (fwd.length()<bigK) return;
		int m=bigK-step;
		if (m<1 || m>31 || step<1)
			throw new IllegalArgumentException("Exact rolling minimizers require 1 <= bigK-step <= 31 and step >= 1");
		long[] canonical=exactCanonicalMers(fwd,m,nc,ncr);
		addExactUnrolledMinimizers(canonical,step,minimizers);
		METRICS.minimizerHashingNanos.add(System.nanoTime()-hashingStartNanos);
	}


	public static void deleteFileOrFolder(File f) throws Exception
	{
		if (f.exists())
		{
			if (f.isDirectory())
			{
				File [] entries = f.listFiles();
				if (entries != null) {for (File entry : entries) {deleteFileOrFolder(entry);}}
			}
			f.delete();
		}
	}
	
		
	public static float[] calculateAvgStd(ConcurrentLinkedQueue<Float> list)
	{
		float[] avgstd = new float[2];
		Iterator iterator = list.iterator();
		while (iterator.hasNext())
		{
			avgstd[0] = avgstd[0] + (float)iterator.next();
		}
		avgstd[0] = avgstd[0] / list.size();
		iterator = list.iterator();
		while (iterator.hasNext())
		{
			float valu = (float)iterator.next();
			avgstd[1] = avgstd[1] + ((avgstd[0] - valu) * (avgstd[0] - valu));
		}
		avgstd[1] = (float)Math.sqrt(avgstd[1] / (list.size()-1));
		return avgstd;
	}
	
	private static String requireValue(String[] args,int index,String option)
	{
		if (index >= args.length){throw new IllegalArgumentException("Missing value after " + option);}
		return args[index];
	}
}

class BloomPopulationPlan
{
	private long singleShardFilters=0L, multiShardFilters=0L;
	private long singleShardInsertions=0L, multiShardInsertions=0L;
	private long totalPhysicalShards=0L;
	private int maximumShardsPerFilter=0;
	private final int bloomWorkers;
	private ExecutorService bloomExecutor;
	private final long[] bloomBatch=new long[BuildCUTTLEFISH.BLOOM_POPULATION_BATCH_LONGS];
	private final byte[] bloomByteBuffer=new byte[BuildCUTTLEFISH.BLOOM_POPULATION_BATCH_LONGS*Long.BYTES];
	private final LongAdder heapInsertionNanos=new LongAdder();
	private final LongAdder mappedParallelInsertionNanos=new LongAdder();
	private final LongAdder heapBatches=new LongAdder();
	private final LongAdder mappedParallelBatches=new LongAdder();
	private long largestPlannedCapacity=0L;
	private int largestPlannedShards=0;
	private final long[] capacityHistogram=new long[7];
	BloomPopulationPlan(int requestedWorkers)
	{
		bloomWorkers=Math.max(1,requestedWorkers);
	}
	private ExecutorService parallelExecutor()
	{
		if (bloomExecutor==null)
		{
			// Create only as many large-stack workers as an active shard batch needs.
			// A fixed pool would eagerly create all configured workers even when a
			// filter has only two or four shards. CallerRuns supplies bounded
			// back-pressure when a future filter has more active shards than workers.
			bloomExecutor=new ThreadPoolExecutor(
				0,bloomWorkers,60L,java.util.concurrent.TimeUnit.SECONDS,
				new SynchronousQueue<Runnable>(),
				new ThreadPoolExecutor.CallerRunsPolicy());
		}
		return bloomExecutor;
	}
	private final IdentityHashMap<GenomesClustersTree,ArrayList<Integer>> capacities =
		new IdentityHashMap<GenomesClustersTree,ArrayList<Integer>>();
	private int filterCount;

	public void addCapacity(GenomesClustersTree node, int capacity)
	{
		ArrayList<Integer> values=capacities.get(node);
		if (values==null)
		{
			values=new ArrayList<Integer>();
			capacities.put(node,values);
		}
		values.add(capacity);
		filterCount++;
	}

	public int getFilterCount() { return filterCount; }
	public long getTotalPhysicalShards() { return totalPhysicalShards; }

	public void populate(GenomesClustersTree root, String outdir, double fpp) throws Exception
	{
		long expectedInsertions=countExpectedInsertions(root,outdir);
		System.out.printf(Locale.ROOT,"\tFilter-local population: %,d expected Bloom insertions%n",expectedInsertions);
		BuildCUTTLEFISH.BloomPopulationProgress progress=
			new BuildCUTTLEFISH.BloomPopulationProgress(expectedInsertions);
		try
		{
			populateNode(root,outdir,fpp,progress);
			progress.finish();
		}
		finally
		{
			if (bloomExecutor!=null)
			{
				bloomExecutor.shutdown();
				if (!bloomExecutor.awaitTermination(1,java.util.concurrent.TimeUnit.MINUTES))
					bloomExecutor.shutdownNow();
			}
		}
		long totalFilters=singleShardFilters+multiShardFilters;
		long totalInsertions=singleShardInsertions+multiShardInsertions;
		double multiFilterPct=totalFilters==0L?0.0:100.0*multiShardFilters/totalFilters;
		double multiInsertionPct=totalInsertions==0L?0.0:100.0*multiShardInsertions/totalInsertions;
		System.out.println("\tBloom shard coverage");
		System.out.printf(Locale.ROOT,"\t  Single-shard filters      = %,d (%.2f%%)%n",singleShardFilters,100.0-multiFilterPct);
		System.out.printf(Locale.ROOT,"\t  Multi-shard filters       = %,d (%.2f%%)%n",multiShardFilters,multiFilterPct);
		System.out.printf(Locale.ROOT,"\t  Single-shard insertions   = %,d (%.2f%%)%n",singleShardInsertions,100.0-multiInsertionPct);
		System.out.printf(Locale.ROOT,"\t  Multi-shard insertions    = %,d (%.2f%%)%n",multiShardInsertions,multiInsertionPct);
		System.out.printf(Locale.ROOT,"\t  Physical shards allocated = %,d; maximum per filter = %,d%n",totalPhysicalShards,maximumShardsPerFilter);
		System.out.printf(Locale.ROOT,"\t  Single-shard backend      = heap-serial with sequential persistence%n");
		System.out.printf(Locale.ROOT,"\t  Multi-shard backend       = mapped-parallel; workers = %,d; super-batch = %,d%n",bloomWorkers,BuildCUTTLEFISH.BLOOM_POPULATION_BATCH_LONGS);
		double heapSeconds=heapInsertionNanos.sum()/1_000_000_000.0;
		double mappedSeconds=mappedParallelInsertionNanos.sum()/1_000_000_000.0;
		System.out.println("\tHybrid Bloom backend profile");
		System.out.printf(Locale.ROOT,"\t  Heap-serial batches       = %,d; insertion time = %.3f s; throughput = %.3f M/s%n",
			heapBatches.sum(),heapSeconds,heapSeconds==0.0?0.0:singleShardInsertions/heapSeconds/1_000_000.0);
		System.out.printf(Locale.ROOT,"\t  Mapped-parallel batches   = %,d; insertion time = %.3f s; throughput = %.3f M/s%n",
			mappedParallelBatches.sum(),mappedSeconds,mappedSeconds==0.0?0.0:multiShardInsertions/mappedSeconds/1_000_000.0);
		System.out.printf(Locale.ROOT,"\t  Largest planned capacity  = %,d elements; %,d physical shards%n",largestPlannedCapacity,largestPlannedShards);
		System.out.printf(Locale.ROOT,"\t  Capacity histogram        = <1M:%,d 1-5M:%,d 5-10M:%,d 10-25M:%,d 25-50M:%,d 50-89.6M:%,d >=89.6M:%,d%n",
			capacityHistogram[0],capacityHistogram[1],capacityHistogram[2],capacityHistogram[3],capacityHistogram[4],capacityHistogram[5],capacityHistogram[6]);
	}

	private long countExpectedInsertions(GenomesClustersTree node, String outdir) throws IOException
	{
		ArrayList<Integer> values=capacities.get(node);
		if (values==null || values.size()!=node.clustersIds.size())
			throw new IllegalStateException("Incomplete Bloom population plan");
		long total=0L;
		for (ArrayList<Integer> genomeIds:node.clustersIds)
			for (int genomeId:genomeIds)
				total=Math.addExact(total,BuildCUTTLEFISH.readLongSetSizeForPlanning(outdir+"genome_"+genomeId));
		for (GenomesClustersTree child:node.childrenList)
			if (child!=null && !child.clustersIds.isEmpty())
				total=Math.addExact(total,countExpectedInsertions(child,outdir));
		return total;
	}

	private void populateNode(GenomesClustersTree node, String outdir, double fpp,
		BuildCUTTLEFISH.BloomPopulationProgress progress) throws Exception
	{
		ArrayList<Integer> values=capacities.get(node);
		if (values==null || values.size()!=node.clustersIds.size())
			throw new IllegalStateException("Incomplete Bloom population plan");
		for (int i=0;i<values.size();i++)
		{
			int capacity=values.get(i);
			long started=System.nanoTime();
			ShardedBloomFilter bf=new ShardedBloomFilter(capacity,fpp,
				outdir+node.clustersBloomIds.get(i));
			int shardCount=bf.getNumShards();
			if ((long)capacity>largestPlannedCapacity)
			{
				largestPlannedCapacity=(long)capacity;
				largestPlannedShards=shardCount;
			}
			int capacityBin=capacity<1_000_000?0:capacity<5_000_000?1:capacity<10_000_000?2:
				capacity<25_000_000?3:capacity<50_000_000?4:capacity<89_617_968?5:6;
			capacityHistogram[capacityBin]++;
			if(shardCount==1)singleShardFilters++;else multiShardFilters++;
			totalPhysicalShards+=shardCount;
			maximumShardsPerFilter=Math.max(maximumShardsPerFilter,shardCount);
			BuildCUTTLEFISH.METRICS.bloomAllocationNanos.add(System.nanoTime()-started);
			try
			{
				ArrayList<Integer> genomeIds=node.clustersIds.get(i);
				for(int genomeId:genomeIds)
				{
					long inserted=BuildCUTTLEFISH.addLongFileToBloom(outdir+"genome_"+genomeId,bf,progress,
						shardCount>1?parallelExecutor():null,bloomBatch,bloomByteBuffer,
						shardCount>1?mappedParallelInsertionNanos:heapInsertionNanos,
						shardCount>1?mappedParallelBatches:heapBatches);
					if(shardCount==1)singleShardInsertions+=inserted;else multiShardInsertions+=inserted;
				}
				if (genomeIds.size()==1) node.clusterElements.set(i,capacity);
				else
				{
					started=System.nanoTime();
					long estimated=bf.getEstimatedDistinctElements();
					BuildCUTTLEFISH.METRICS.bloomEstimationNanos.add(System.nanoTime()-started);
					node.clusterElements.set(i,
						Math.round(((float)capacity+(float)estimated)/2));
				}
			}
			finally
			{
				started=System.nanoTime();
				bf.close();
				BuildCUTTLEFISH.METRICS.bloomCloseNanos.add(System.nanoTime()-started);
			}
		}
		for (GenomesClustersTree child:node.childrenList)
			if (child!=null && !child.clustersIds.isEmpty())
				populateNode(child,outdir,fpp,progress);
	}
}

class GenomesClustersTree implements Serializable
{
	private static final AtomicLong BLOOM_ID = new AtomicLong();
	private static final ThreadLocal<Integer> CONSTRUCTION_DEPTH =
		ThreadLocal.withInitial(() -> Integer.valueOf(0));
	private static final AtomicInteger MAX_CONSTRUCTION_DEPTH = new AtomicInteger(0);
	public static void resetConstructionDepthDiagnostics()
	{
		CONSTRUCTION_DEPTH.set(Integer.valueOf(0));
		MAX_CONSTRUCTION_DEPTH.set(0);
	}
	public static int getMaximumConstructionDepth()
	{
		return MAX_CONSTRUCTION_DEPTH.get();
	}
	public ArrayList<ArrayList<Integer>> clustersIds;
	public ArrayList<String> clustersBloomIds;
	public ArrayList<GenomesClustersTree> childrenList;
	public ArrayList<Integer> clusterElements;

	public GenomesClustersTree()
	{
		clustersIds = new ArrayList<ArrayList<Integer>>();
		clustersBloomIds = new ArrayList<String>();
		childrenList = new ArrayList<GenomesClustersTree>();
		clusterElements = new ArrayList<Integer>();
	}
	
	public GenomesClustersTree(ArrayList<Integer> genomesIds, ConcurrentHashMap<Long,Double> similarityMatrixMap, ArrayList<BuildCUTTLEFISH.GenomeSketch> genomesHLL, int numThreads, String outdir, AtomicInteger processedGenomes, BufferedWriter bw, AtomicLong numNodes, double fpp, int maxKmersPerSpecies, int kappa, double globalAvgSimilarity, double globalStdSimilarity) throws Exception
	{
		this(genomesIds, similarityMatrixMap, genomesHLL, numThreads, outdir, processedGenomes, bw, numNodes, fpp, maxKmersPerSpecies, kappa, globalAvgSimilarity, globalStdSimilarity, null);
	}
	
	public GenomesClustersTree(ArrayList<Integer> genomesIds, ConcurrentHashMap<Long,Double> similarityMatrixMap, ArrayList<BuildCUTTLEFISH.GenomeSketch> genomesHLL, int numThreads, String outdir, AtomicInteger processedGenomes, BufferedWriter bw, AtomicLong numNodes, double fpp, int maxKmersPerSpecies, int kappa, double globalAvgSimilarity, double globalStdSimilarity, BloomPopulationPlan bloomPlan) throws Exception
	{	
		int constructionDepth=CONSTRUCTION_DEPTH.get().intValue()+1;
		CONSTRUCTION_DEPTH.set(Integer.valueOf(constructionDepth));
		MAX_CONSTRUCTION_DEPTH.accumulateAndGet(constructionDepth,Math::max);
		clustersIds = new ArrayList<ArrayList<Integer>>();
		clustersBloomIds = new ArrayList<String>();
		childrenList = new ArrayList<GenomesClustersTree>();
		clusterElements = new ArrayList<Integer>();
		long clusteringStartNanos=System.nanoTime();
		ArrayList<ArrayList<Integer>> clustersA = BuildCUTTLEFISH.clusterGenomes(genomesIds, similarityMatrixMap, genomesHLL, numThreads, kappa, globalAvgSimilarity, globalStdSimilarity);
		BuildCUTTLEFISH.METRICS.clusteringNanos.add(System.nanoTime()-clusteringStartNanos);
		if (clustersA.size()==1)
		{
			for (int p=0; p<clustersA.get(0).size(); p++) {numNodes.incrementAndGet();} //System.out.println("\tNode no. "+numNodes.get()+" singleton of size "+clustersA.get(0).size());
			if (genomesIds.size()==1) {bw.write(""+genomesIds.get(0));}
			else {bw.write("("); for (int j=0; j<genomesIds.size()-1; j++) {bw.write(genomesIds.get(j)+",");} bw.write(genomesIds.get(genomesIds.size()-1)+")");}
			ArrayList<Integer> cluster = clustersA.get(0);
			for (int j=0; j<cluster.size(); j++)
			{
				ArrayList<Integer> genomeName = new ArrayList<Integer>();
				int genomeId = cluster.get(j);
				genomeName.add(genomeId);
				processedGenomes.incrementAndGet();
				clustersIds.add(genomeName);
				String bloomName = "bloom_" + BLOOM_ID.getAndIncrement() + ".blo";
				long genomeCardinalityLong = BuildCUTTLEFISH.readLongSetSize(outdir+"genome_"+genomeId);
				if (genomeCardinalityLong > Integer.MAX_VALUE)
					throw new IOException("Genome minimizer count exceeds supported range: " + genomeCardinalityLong);
				int genomeCardinality = (int)genomeCardinalityLong;
				if (bloomPlan != null)
				{
					bloomPlan.addCapacity(this, genomeCardinality);
					clusterElements.add(genomeCardinality);
				}
				else
				{
					long bloomAllocationStartNanos=System.nanoTime();
					ShardedBloomFilter bf = new ShardedBloomFilter(genomeCardinality,fpp,outdir+bloomName);
					BuildCUTTLEFISH.METRICS.bloomAllocationNanos.add(System.nanoTime()-bloomAllocationStartNanos);
					BuildCUTTLEFISH.addLongFileToBloom(outdir+"genome_"+genomeId, bf);
					clusterElements.add(genomeCardinality);
					long bloomCloseStartNanos=System.nanoTime();
					bf.close();
					BuildCUTTLEFISH.METRICS.bloomCloseNanos.add(System.nanoTime()-bloomCloseStartNanos);
				}
				clustersBloomIds.add(bloomName);
				childrenList.add(new GenomesClustersTree());
				if (processedGenomes.get()%1000==0)
				{
					System.out.printf(Locale.ROOT,
						"	Tree progress: %,d/%,d genomes assigned; %,d nodes created%n",
						processedGenomes.get(),genomesHLL.size(),numNodes.get());
				}
			}
		}
		else
		{
			bw.write("(");
			LinkedList<ArrayList<Integer>> clusters = new LinkedList<ArrayList<Integer>>();
    			clusters.addAll(clustersA);
			while(clusters.size()>0)
			{
				numNodes.incrementAndGet(); //System.out.println("\tNode no. "+numNodes.get());
				ArrayList<Integer> cluster = clusters.removeFirst();
				int maxElements = 0;
				long sketchCapacityStartNanos=System.nanoTime();
				long clusterSize=genomesHLL.get(cluster.get(0)).cardinality();
				if(clusterSize>maxKmersPerSpecies)clusterSize=maxKmersPerSpecies;
				Union capacityUnion=SetOperation.builder().setNominalEntries(BuildCUTTLEFISH.THETA_K).setSeed(BuildCUTTLEFISH.THETA_SEED).buildUnion();
				for(int j=0;j<cluster.size();j++){
					capacityUnion.union(genomesHLL.get(cluster.get(j)).theta);
					long cs=Math.max(0L,Math.round(capacityUnion.getResult().getEstimate()));
					if(cs<=maxKmersPerSpecies){maxElements=j;clusterSize=cs;}else break;
				}
				BuildCUTTLEFISH.METRICS.sketchCapacityNanos.add(System.nanoTime()-sketchCapacityStartNanos);
				//clusterSize=(int)(clusterSize+0.02d*clusterSize); if (clusterSize>Integer.MAX_VALUE) {clusterSize=Integer.MAX_VALUE;}
				//System.out.print("writing bloom.. ");
				String bloomName = "bloom_" + BLOOM_ID.getAndIncrement() + ".blo";
		//bloomName=stron;
				ShardedBloomFilter bf = null;
				if (bloomPlan != null) bloomPlan.addCapacity(this, (int)clusterSize);
				else
				{
					long bloomAllocationStartNanos=System.nanoTime();
					bf = new ShardedBloomFilter((int)clusterSize,fpp,outdir+bloomName);
					BuildCUTTLEFISH.METRICS.bloomAllocationNanos.add(System.nanoTime()-bloomAllocationStartNanos);
				}
				ArrayList<Integer> splitCluster = new ArrayList<Integer>();
				for (int j=maxElements+1; j<cluster.size(); j++)
				{
					splitCluster.add(cluster.get(j));
				}
				if (splitCluster.size()>0) {clusters.add(splitCluster);}
				ArrayList<Integer> clusterIds = new ArrayList<Integer>();
				for (int j=0; j<=maxElements; j++)
				{
					int genomeId = cluster.get(j);
					clusterIds.add(genomeId);
					if (bloomPlan == null) BuildCUTTLEFISH.addLongFileToBloom(outdir+"genome_"+genomeId, bf);
				}
				clustersIds.add(clusterIds);
				if (bloomPlan != null) clusterElements.add((int)clusterSize);
				else
				{
					long bloomEstimationStartNanos=System.nanoTime();
					long estimatedDistinctElements=bf.getEstimatedDistinctElements();
					BuildCUTTLEFISH.METRICS.bloomEstimationNanos.add(System.nanoTime()-bloomEstimationStartNanos);
					clusterElements.add(Math.round(((float)clusterSize + (float)estimatedDistinctElements) / 2));
					long bloomCloseStartNanos=System.nanoTime();
					bf.close();
					BuildCUTTLEFISH.METRICS.bloomCloseNanos.add(System.nanoTime()-bloomCloseStartNanos);
				}
				clustersBloomIds.add(bloomName); //System.out.println(" done");
				if (clusterIds.size()>1) childrenList.add(new GenomesClustersTree(clusterIds, similarityMatrixMap, genomesHLL, numThreads, outdir, processedGenomes, bw, numNodes, fpp, maxKmersPerSpecies, kappa, globalAvgSimilarity, globalStdSimilarity, bloomPlan));
				else 
				{ //System.out.println("\t singleton node ");
					childrenList.add(new GenomesClustersTree());
					bw.write(""+clusterIds.get(0));
					processedGenomes.incrementAndGet();
					if (processedGenomes.get()%1000==0)
					{
						System.out.printf(Locale.ROOT,
							"	Tree progress: %,d/%,d genomes assigned; %,d nodes created%n",
							processedGenomes.get(),genomesHLL.size(),numNodes.get());
					}
				}
		//if (clusterIds.size()>1) bw.write(stron);
				if (clusters.size()>0) bw.write(",");
			}
			bw.write(")"); bw.newLine();
		}
		CONSTRUCTION_DEPTH.set(Integer.valueOf(constructionDepth-1));
	}
	
}


