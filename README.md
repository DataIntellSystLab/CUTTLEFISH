CUTTLEFISH
==========

Classification by U-mer Tuples of Taxa with Likelihood Estimated Filters Incrementally Structured under a Hierarchy

OVERVIEW
--------

CUTTLEFISH is a portable, multiplatform metagenomics classifier designed for offline analysis of nanopore reads on systems with constrained memory, including laptops, tablets, smartphones, and other edge-computing devices. It represents reference genomes through distinct canonical k-mer minimizers, organizes related taxa in a similarity-derived tree, and stores each branch in independently sized, optionally sharded Bloom filters.

CUTTLEFISH supports a configurable balance among classification accuracy, index size, construction time, runtime, and memory use. During classification, Bloom filters that fit within a user-defined RAM budget are loaded into memory, while the remaining filters are accessed through read-only memory mapping. Consequently, an index may be used even when it is larger than the available Java heap, subject to sufficient local storage.

The implementation has been tested with viral, mitochondrial, fungal, bacterial, and antimicrobial-resistance reference collections. On held-out real genomes with nanopore-like sequencing reads, median classification performance was comparable to Kraken2 after selecting a confidence threshold appropriate for the Bloom-filter false-positive probability, while variability across species was substantially lower.

The reference databases' indices are available at [archive URL/DOI].

IMPLEMENTATION
--------------

CUTTLEFISH is implemented in Java and uses multithreading for database construction and read classification. The current implementation targets Java 26 and uses the Java standard library together with:

* Apache Commons Math
* Apache DataSketches
* fastutil
* hash4j
* stream-lib

Reference genomes are represented as sets of distinct canonical minimizers. By default, CUTTLEFISH uses k = 29. For each k-base window, it evaluates adjacent canonical m-mers, where m = k - step, and retains the minimum. The default step is floor(sqrt(k)) rounded down to the nearest even integer. Valid k values are odd integers from 11 through 37. With the currently supported parameterization, step is 2, 4, or 6 and m does not exceed 31, allowing collision-free two-bit encoding in a Java long.

Genome similarity is estimated using Theta sketches and SuperMinHash signatures. A deterministic, two-pass multi-leader clustering procedure recursively organizes genomes into a taxonomic hierarchy. Internal and terminal branches are represented by Bloom filters. The Bloom-filter false-positive probability is configurable when building an index and controls a principal storage-versus-accuracy trade-off. The default is 1e-5; a larger value, such as 0.0075, produces a smaller index and usually a faster build, but may require a more conservative classification-confidence threshold.

During classification, the program loads complete Bloom filters into the Java heap until the requested filter-memory budget is exhausted. Remaining filters are opened through read-only memory mapping. Read admission is also bounded according to thread count, estimated read size, and heap pressure. Species-level distinct minimizers are counted using hash4j UltraLogLog when the estimated counter bank fits the available-memory criterion; otherwise, Apache DataSketches HLL_4 counters are created lazily.

FILES
-----

The source distribution is expected to contain at least:

* BuildCUTTLEFISH.java     database-index builder and hierarchy construction
* CUTTLEFISH.java          FASTQ read classifier and summary writer
* ShardedBloomFilter.java  sharded Bloom-filter storage, hashing, RAM loading, and memory mapping
* cuttlefish_jars/         required third-party JAR files

A completed database directory contains metadata.txt, taxonomicClusterTree.ser, info.txt, and the Bloom-filter files. CUTTLEFISH validates database-format and minimizer metadata before classification and rejects incomplete or incompatible indices.

COMPILATION
-----------

Linux or macOS:

    javac -cp ".:cuttlefish_jars/*" BuildCUTTLEFISH.java CUTTLEFISH.java ShardedBloomFilter.java

Windows Command Prompt or PowerShell:

    javac -cp ".;cuttlefish_jars/*" BuildCUTTLEFISH.java CUTTLEFISH.java ShardedBloomFilter.java

The classpath separator is ':' on Linux/macOS and ';' on Windows.

To display the command-line help supplied by either program:

    java -cp ".:cuttlefish_jars/*" BuildCUTTLEFISH --help
    java -cp ".:cuttlefish_jars/*" CUTTLEFISH --help

Replace ':' with ';' on Windows.

BUILDING A DATABASE INDEX
-------------------------

Basic form:

    java [JVM options] -cp ".:cuttlefish_jars/*" BuildCUTTLEFISH \
        -f references.fasta.gz [builder options]

Example to build a large reference collection (e.g., 10,000 bacterial genomes) with loose false positive rate:

    java -Xmx32G -Xss128m -cp ".:cuttlefish_jars/*" BuildCUTTLEFISH \
        -f bacterial_genomes.fna.gz --e 0.01

Example to build a small reference collection (e.g., 10,000 viral genomes) with strict false positive rate:

    java -Xmx16G -Xss128m -cp ".:cuttlefish_jars/*" BuildCUTTLEFISH \
        -f viral_genomes.fna.gz -e 0.00001

Builder arguments:

    -f PATH       Required. Input may be a FASTA file, gzip-compressed FASTA, or
                  a directory containing genome FASTA files. Directory input is
                  searched recursively for .fna, .fa, .fasta, and .fas files,
                  optionally ending in .gz. Protein, translated-CDS, RNA, and
                  CDS-derived files are excluded by filename.

    -k INTEGER    k-mer window length.
                  Default: 29. Accepted range: 11 through 37.
                  If an even value is supplied, the program increments it to the
                  next odd value before enforcing the range. The minimizer m-mer
                  length is m = k - step, where step is floor(sqrt(k)) rounded
                  down to an even integer.

    -m INTEGER    Process only the first INTEGER genomes.
                  Default: all genomes. Values less than or equal to zero revert
                  to all genomes.

    -t INTEGER    Number of worker threads.
                  Default: available processors minus one, with a minimum of one.
                  Nonpositive values revert to the default.

    -e DECIMAL    Target Bloom-filter false-positive probability.
                  Default: 0.00001 (1e-5). Valid values are strictly between 0
                  and 1. Invalid values revert to 0.00001. Larger values reduce
                  Bloom-filter and overall index size, but generally require a
                  more conservative classification threshold.

    --clustering-seed LONG
                  Deterministic clustering seed. Default: 9001. Reusing the same
                  input, parameters, software version, and seed makes the seeded
                  clustering decisions reproducible.

    -h, --help    Print command-line help and exit.

The builder automatically names the output directory from the input basename as:

    INPUT_BASENAME_CUTTLEFISHdb/

Warning: if that directory already exists, the builder deletes it before starting.
An initial metadata.txt is written with buildComplete=false and is changed to true
only after successful construction and serialization. The completed metadata also
records k, step, false-positive probability, clustering seed, genome count, Bloom-
filter count, and physical-shard count. The distinct-counting precision option
(-l) belongs to classification summaries and is not a database-build parameter.

Practical guidance:

* Bloom filters are sized independently. A logical filter that exceeds one physical filter's supported capacity is partitioned into a power-of-two number of shards, up to 255. Each shard uses a prime number of bits and MurmurHash2 double hashing; a seed-0 hash selects the shard. Single-shard filters are built in heap and persisted sequentially, while multi-shard filters use mapped, shard-parallel insertion.

* Use -Xmx to set the maximum Java heap appropriate for the host. Index construction is designed to operate with bounded memory, but large collections need adequate temporary storage and may benefit from a larger heap.
* Use -t conservatively on mobile or thermally constrained hardware. More threads can increase throughput but also increase transient memory use and I/O pressure.
* Keep the default Bloom-filter false-positive probability of 1e-5 when classification consistency is the priority.
* Increase the Bloom-filter false-positive probability when a smaller index is required. For example, 0.0075 substantially reduces index size, but the classification threshold should then be tuned or chosen conservatively.
* Database construction can be slower than classification-oriented server tools, but indices are built once and can be distributed as precomputed resources.
* Do not modify or remove metadata.txt. A database is accepted only when buildComplete=true and all format identifiers match the classifier.

CLASSIFYING READS
-----------------

Basic form:

    java [JVM options] -cp ".:cuttlefish_jars/*" CUTTLEFISH \
        -d DATABASE_DIRECTORY -f READS.fastq.gz [options]

Example using a 2 GB Java heap and retaining every positive match:

    java -Xmx2G -Xss128m -cp ".:cuttlefish_jars/*" CUTTLEFISH \
        -d bacterial_genomes_CUTTLEFISHdb \
        -f reads.fastq.gz \
        -p 0

Required classifier arguments:

    -d PATH       CUTTLEFISH database-index directory.
    -f PATH       Input FASTQ file. Plain-text and gzip-compressed FASTQ are supported.

Classifier options:

    -m NUMBER     Classify exactly the first NUMBER reads.
                  Default: all reads.
                  A negative value is treated as no limit.

    -t INTEGER    Number of worker threads.
                  Default: available processors minus one, with a minimum of one.
                  Nonpositive values revert to the default.

    -p DECIMAL    Acceptance threshold.
                  Default: 0.5.
                  For values below 1, this is the terminal confidence threshold,
                  where confidence = score / (1 + score).
                  For values greater than or equal to 1, the value is interpreted
                  as a minimum number of matching read minimizers.
                  A negative value is rejected and reset to 0.5.
                  Use 0 to retain every classification with at least one positive
                  database match for downstream threshold analysis.

    -b DECIMAL    Fraction of the JVM's maximum heap reserved for loading complete
                  Bloom filters into RAM.
                  Default: 0.5.
                  Valid range: 0.05 through 0.95, inclusive.
                  Invalid values revert to 0.5.
                  Filters are loaded in RAM while this budget permits; all remaining
                  filters are accessed through read-only memory mapping.

    -l INTEGER    Log2 precision for species-level distinct-minimizer counters.
                  Default: 12.
                  Valid range: 10 through 31, inclusive.
                  Invalid values revert to 12.
                  Greater precision increases counter memory use. CUTTLEFISH selects
                  UltraLogLog or lazy HLL_4 storage automatically according to the
                  estimated memory requirement.

    -h, --help    Print command-line help and exit.

Memory examples:

1. Favor low heap usage and disk-backed access:

    java -Xmx512M -cp ".:cuttlefish_jars/*" CUTTLEFISH \
        -d DATABASE_DIRECTORY -f READS.fastq.gz -b 0.25 -t 2

2. Load more filters into RAM when memory is available:

    java -Xmx16G -cp ".:cuttlefish_jars/*" CUTTLEFISH \
        -d DATABASE_DIRECTORY -f READS.fastq.gz -b 0.70 -t 15

A larger -b value is not always faster. Leaving too little heap and operating-system headroom can increase Java heap-pressure backpressure and reduce the page cache available to memory-mapped filters. Tune -Xmx, -b, and -t together for the target device and storage system.

Threshold guidance:

* The default classification threshold is 0.5.
* The appropriate threshold depends partly on the Bloom-filter false-positive probability used to construct the index.
* For an index built at 1e-5, lower positive thresholds may improve sensitivity while retaining high specificity.
* A compact index built with a larger false-positive probability generally requires a higher threshold to suppress Bloom-filter false positives.
* For a new reference collection or use case, run with -p 0, retain all positive matches, and select an operating threshold using representative validation data.

OUTPUT
------

For an input named READS.fastq.gz, CUTTLEFISH writes two CSV files beside the input:

    READS.fastq.gz_CUTTLEFISH_classReads.csv
    READS.fastq.gz_CUTTLEFISH_classGenos.csv

Per-read output columns:

    readId,genoId,genoMatches|readMinimizers,percMatch,confidenceScore

* readId: FASTQ read header.
* genoId: assigned genome/species name; '?' indicates rejection or no positive match.
* genoMatches|readMinimizers: matching minimizers and total distinct query minimizers.
* percMatch: fraction of query minimizers matching the selected terminal genome.
* confidenceScore: bounded score calculated as S/(1+S).

Per-genome output columns:

    genoId,genoName,genoCoverage,kmerMinimizDepth

* genoCoverage: estimated percentage of the reference genome's distinct minimizers observed in accepted reads, capped at 100%.
* kmerMinimizDepth: accepted minimizer-hit depth divided by the estimated number of observed distinct minimizers.

INPUT NOTES
-----------

* FASTA reference files and FASTQ read files may be gzip-compressed. A multi-record FASTA treats each FASTA record as a reference entry; with directory input, each accepted FASTA file is treated as a genome and its contigs are combined.
* FASTQ parsing accepts wrapped sequence and quality lines and attempts lenient resynchronization after malformed records.
* Characters other than A, C, G, T, or U reset minimizer generation, preventing minimizers from spanning ambiguous or invalid positions.
* Database indices are tied to their recorded minimizer and Bloom-filter formats. Rebuild an index if the classifier reports incompatible metadata.

REPRODUCIBILITY AND PORTABILITY
-------------------------------

Clustering is reproducible for a fixed seed. The software is Java-based to facilitate deployment across operating systems and hardware architectures. Actual runtime depends on read length, thread count, index size, Java heap, the RAM-loading fraction, local-storage latency, filesystem contention, and operating-system page-cache behavior. Fully RAM-resident indices are generally faster and more stable; larger memory-mapped indices trade speed for the ability to operate under strict memory limits.

LICENSE
-------

Apache License 2.0.

CITATION
--------

If you use CUTTLEFISH, please cite:

Marini S, Boucher C, Ruiz J, Prosperi M. CUTTLEFISH, a portable metagenomics classifier based on Bloom filter hierarchies. [Journal, year, DOI to be added].
