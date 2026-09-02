//Classification by U-mer Tuples of Taxa with Likelihood Estimated Filters Incrementally Structured under a Hierarchy (CUTTLEFISH)

/* 
del *.class
javac -cp ".;cuttlefish_jars/*" BuildCUTTLEFISH.java
javac -cp ".;cuttlefish_jars/*" CUTTLEFISH.java
java -Xmx23G -Xss128m -cp ".;cuttlefish_jars/*" BuildCUTTLEFISH -f bacteria_who_2023.fasta.gz -t 23 -m 125
java -Xmx23G -Xss128m -cp ".;cuttlefish_jars/*" BuildCUTTLEFISH -f viruses_refseq_2023.fna.gz -t 23 -m 2222
java -Xmx2G -Xss128m -cp ".;cuttlefish_jars/*" CUTTLEFISH -d bacteria_who_2023_CUTTLEFISHdb -f bacterial_sampled.fastq.gz -p 0 -t 23
java -Xmx2G -Xss128m -cp ".;cuttlefish_jars/*" CUTTLEFISH -d viruses_refseq_2023_CUTTLEFISHdb -f simulmix.fastq -p 0 -t 23
*/

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
import java.util.zip.*;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import com.dynatrace.hash4j.distinctcount.UltraLogLog;
import org.apache.datasketches.hll.HllSketch;
import org.apache.datasketches.hll.TgtHllType;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.longs.LongIterator;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;

public class CUTTLEFISH
{
	public static void main(String[] args) throws Exception
	{
		long timeZero = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long elapsedTime = System.currentTimeMillis() - startTime;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		System.out.println();
		System.out.println("Max RAM allocated is "+allram/1048576+" MB");
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		final int DEFAULT_BUFFER_SIZE = 16384;
		int cores = Runtime.getRuntime().availableProcessors()-1;
		long max_reads = Long.MAX_VALUE;
		float desired_confidence_or_minFreq = 0.5f;
		double rambudg=0.5d;
		int hll_log2m = 12;
		String dbdir1="";
		String fastqfile="";
		if (args.length<=0) {System.out.println("type java CUTTLEFISH h [or help, -h, -help] for instructions"); System.exit(0);}
		for (int t=0; t<args.length; t++)
		{
			switch (args[t])
			{
				case "-h":
				case "--help":
					System.out.println(
						"Usage:\r\n" +
						"\tjava CUTTLEFISH -f <FASTQ file> -d <database folder> [options]\r\n\r\n" +
						"Required:\r\n" +
						"\t-f <path>       FASTQ file, optionally gzipped\r\n" +
						"\t-d <path>       CUTTLEFISH database folder\r\n\r\n" +
						"Options:\r\n" +
						"\t-m <number>     classify exactly the first number of reads\r\n" +
						"\t-t <integer>    number of worker threads, default is all processors minus one\r\n" +
						"\t-p <decimal>    confidence or minimum-frequency threshold, default is 0.5 confidence\r\n" +
						"\t-b <decimal>    portion of RAM out of total allocated for the JVM for loading Bloom the filters (with 0.05 tolerance), default is 0.5\r\n" +
						"\t-l <integer>    distinct-counter precision, default is 12\r\n" +
						"\t-h, --help      print this help"
					);
					return;

				case "-m":
					max_reads = Long.parseLong(requireValue(args, ++t, "-m"));
					if (max_reads < 0)
					{
						max_reads = Long.MAX_VALUE;
					}
					break;

				case "-p":
					desired_confidence_or_minFreq = Float.parseFloat(
						requireValue(args, ++t, "-p")
					);
					if (desired_confidence_or_minFreq < 0)
					{
						desired_confidence_or_minFreq = 0.5f;
						System.out.println(
							"Invalid confidence or frequency threshold; using 0.5."
						);
					}
					break;
					
				case "-b":
				rambudg = Float.parseFloat(
					requireValue(args, ++t, "-b")
				);
				if (rambudg<0.05 || rambudg>0.95)
				{
					rambudg = 0.5f;
					System.out.println(
						"Invalid % RAM budget for loading Bloom filters; using 0.5."
					);
				}
				break;

				case "-l":
					hll_log2m = Integer.parseInt(
						requireValue(args, ++t, "-l")
					);
					if (hll_log2m < 10 || hll_log2m > 31)
					{
						hll_log2m = 12;
						System.out.println(
							"Invalid log2m value; using 12."
						);
					}
					break;

				case "-t":
					cores = Integer.parseInt(
						requireValue(args, ++t, "-t")
					);
					if (cores <= 0)
					{
						cores = Runtime.getRuntime().availableProcessors()-1;
						System.out.println(
							"Invalid thread count; using all processors minus one."
						);
					}
					break;

				case "-d":
					dbdir1 = requireValue(args, ++t, "-d");
					break;

				case "-f":
					fastqfile = requireValue(args, ++t, "-f");
					break;

				default:
					throw new IllegalArgumentException(
						"Unknown option: " + args[t] + ". Use -h for help."
					);
			}
		}
		if (cores<=0) cores=1;
		if (dbdir1.equals("")) {System.out.println("Please specify a CUTTLEFISH database with -d database_index_folder"); System.exit(0);}
		final String dbdir = dbdir1;
		if (fastqfile.equals("")) {System.out.println("Please specify a FASTQ input file with -f fastq_file (can be gzipped)"); System.exit(0);}
		Properties p = loadAndValidateMetadata(dbdir);
		System.out.println("Running program on "+cores+" threads");
		
		System.out.println("Reading taxonomic clusters' tree and ancillary files");
		startTime = System.currentTimeMillis();
		InputStream is = new FileInputStream(new File(dbdir+"/taxonomicClusterTree.ser"));
		ObjectInputStream ois = new ObjectInputStream(is);
		GenomesClustersTree gct = (GenomesClustersTree)ois.readObject();
		ois.close(); is.close();
		HashMap<Integer,CuttlefishSpecies> panGenomeSummaries = new HashMap<Integer,CuttlefishSpecies>();
		BufferedReader br = new BufferedReader(new FileReader(dbdir+"/info.txt"));
		String line = br.readLine();
		line = br.readLine();
		ArrayList<Integer> genoIdSet = new ArrayList<Integer>();
		while(line!=null)
		{
			String[] item = line.split(",");
			CuttlefishSpecies cs = new CuttlefishSpecies();
			cs.totKmerMinimiz = Long.parseLong(item[1]);
			cs.speciesName = item[2];
			int genoId = Integer.parseInt(item[0]);
			panGenomeSummaries.put(genoId,cs);
			genoIdSet.add(genoId);
			line = br.readLine();
		}
		br.close();
		System.out.println("\tThe database contains "+ panGenomeSummaries.size() +" species");
		long estimatedUltraBytes=(long)panGenomeSummaries.size()*((1L<<hll_log2m)+64L);
		long fastCounterLimit=Runtime.getRuntime().maxMemory()/8L;
		final boolean useUltraCounters=estimatedUltraBytes<=fastCounterLimit;
		CounterFactory.configure(useUltraCounters,hll_log2m);
		System.out.printf(Locale.ROOT,
			"\tDistinct-minimizer counters: AUTO selected %s, precision %d%s%n",
			useUltraCounters?"hash4j UltraLogLog":"lazy DataSketches HLL_4",hll_log2m,
			useUltraCounters?String.format(Locale.ROOT,"; estimated bank %.1f MiB",estimatedUltraBytes/1048576.0):"");
		System.out.println("\tDistinct-counter input mixer: seeded Murmur3 fmix64");
		System.out.println("\tRead-classification traversal: primitive minimizer array");
		for (int genoId:genoIdSet)
		{
			CuttlefishSpecies cs=panGenomeSummaries.get(genoId);
			cs.curKmerMinimiz=useUltraCounters?CounterFactory.create():null;
		}


		int k = parsePositiveMetadataInt(p,"k");
		int stp = parsePositiveMetadataInt(p,"step");
		if (stp>=k) throw new IOException("Invalid database metadata: step must be smaller than k");
		final int step = stp;
		System.out.println("Length of kmers is "+k+" (minimizers are k-"+step+", i.e. "+(k-step)+")");
		
		HashMap<String,ShardedBloomFilter> bfh = new HashMap<String,ShardedBloomFilter>();
		BufferedReader fastqReader = null;
		BufferedWriter bw = null;
		BufferedWriter pw = null;
		ExecutorService executor = null;
		Throwable primaryFailure = null;
		try
		{
		long memoryBudgetBytes = Math.round(rambudg*(double)Runtime.getRuntime().maxMemory());
		System.out.println("\tOpening Bloom filters (up to "+rambudg+" of available RAM)");
		AtomicInteger filtersLoadedInRAM = new AtomicInteger(0);
		AtomicInteger filtersLoadedMMAP = new AtomicInteger(0);
		long remainingBudgetBytes = openFilters(bfh,gct,dbdir,memoryBudgetBytes,filtersLoadedInRAM,filtersLoadedMMAP);
		System.out.println("\t\t"+bfh.size()+" filters opened");
		System.out.println("\t\t"+filtersLoadedInRAM.get()+" RAM-resident filters");
		System.out.println("\t\t"+filtersLoadedMMAP.get()+" memory-mapped filters");
		System.out.printf(Locale.ROOT,"\t\t%.1f MiB RAM budget remaining%n",remainingBudgetBytes/1048576.0);
		System.gc();
		elapsedTime = System.currentTimeMillis() - startTime;
		System.out.print("\tCore and support files read in "+elapsedTime/1000+" s");
		usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		System.out.println("; used RAM "+usedram/1048576+" MB");
		
		System.out.println("Processing FASTQ file");
		startTime = System.currentTimeMillis();
		final long [][] pp = BuildCUTTLEFISH.precomputedPowsForNucs(4,31);
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
		fastqReader = new BufferedReader(new FileReader(fastqfile),DEFAULT_BUFFER_SIZE);
		if(fastqfile.endsWith(".gz")) fastqReader=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fastqfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		ArrayDeque<String> pendingFastqLines = new ArrayDeque<String>();
		line = readFastqLine(fastqReader, pendingFastqLines);
		bw = new BufferedWriter(new FileWriter(fastqfile+"_CUTTLEFISH_classReads.csv"), DEFAULT_BUFFER_SIZE);
		pw = new BufferedWriter(new FileWriter(fastqfile+"_CUTTLEFISH_classGenos.csv"), DEFAULT_BUFFER_SIZE);
		bw.write("readId,genoId,genoMatches|readMinimizers,percMatch,confidenceScore"); bw.newLine();
		pw.write("genoId,genoName,genoCoverage,kmerMinimizDepth"); pw.newLine();
		executor = Executors.newFixedThreadPool(cores);
		CompletionService<CuttlefishResult> completionService =
		new ExecutorCompletionService<CuttlefishResult>(executor);
		final int MAX_TASKS_IN_FLIGHT = Math.max(cores * 2, 8);
		final long PROCESSING_BUDGET_BYTES=Math.max(64L<<20,
			Math.min(1L<<30,Runtime.getRuntime().maxMemory()*15L/100L));
		long estimatedBytesInFlight=0L;
		long peakEstimatedBytesInFlight=0L;
		long maximumReadLength=0L;
		long totalReadLength=0L;
		long taskWindowWaits=0L;
		long processingBudgetWaits=0L;
		long heapPressureWaits=0L;
		int peakTasksInFlight=0;
		int tasks_in_flight = 0;
		System.out.printf(Locale.ROOT,"\tAdaptive read admission: %.1f MiB budget; at most %,d reads in flight%n",
			PROCESSING_BUDGET_BYTES/1048576.0,MAX_TASKS_IN_FLIGHT);
		long reads_processed = 0L;
		long reads_submitted = 0L;
		//ConcurrentLinkedQueue<CuttlefishResult> processed = new ConcurrentLinkedQueue<CuttlefishResult>();
		while(line!=null && reads_submitted<max_reads)
		{
			boolean nextHeaderAlreadyRead = false;
			if (!line.startsWith("@"))
			{
				// Lenient resynchronization: skip stray lines until a plausible
				// header is found instead of terminating the entire run.
				System.err.println("\tWarning: skipping non-header FASTQ line while resynchronizing: " + line.substring(0, Math.min(line.length(), 10)));
				do
				{
					line = readFastqLine(fastqReader, pendingFastqLines);
				}
				while (line != null && !line.startsWith("@"));
				if (line == null) break;
			}
			String read_header = line;
			StringBuilder seqBuilder = new StringBuilder();
			String separator = null;
			while ((line = readFastqLine(fastqReader, pendingFastqLines)) != null)
			{
				if (line.startsWith("+"))
				{
					separator = line;
					break;
				}
				seqBuilder.append(line.trim());
			}
			if (separator == null)
			{
				System.err.println("\tWarning: incomplete FASTQ record for " + read_header.substring(0, Math.min(read_header.length(), 10)) +
					": missing '+' separator; stopping at end of file.");
				break;
			}
			String read_sequence = seqBuilder.toString();
			StringBuilder qualBuilder = new StringBuilder();
			while (qualBuilder.length() < read_sequence.length())
			{
				String qline = readFastqLine(fastqReader, pendingFastqLines);
				if (qline == null)
				{
					break;
				}

				// '@' is legal quality data. Treat it as a header only when the
				// following buffered lines look like a complete FASTQ sequence
				// section ending in '+'. This permits wrapped quality strings that
				// themselves contain lines beginning with '@'.
				if (qline.startsWith("@") &&
					isPlausibleNextFastqHeader(fastqReader, pendingFastqLines))
				{
					line = qline;
					nextHeaderAlreadyRead = true;
					break;
				}
				qualBuilder.append(qline);

				// Once quality is longer than sequence, boundaries are inherently
				// ambiguous. Preserve and process this read, warn below, and resume
				// from the next physical line (the recoverable single-line case).
				if (qualBuilder.length() > read_sequence.length())
				{
					break;
				}
			}
			String quality = qualBuilder.toString();
			if (quality.length() != read_sequence.length())
			{
				System.err.println("\tWarning: sequence and quality lengths differ for " + read_header.substring(0, Math.min(read_header.length(), 10)) +
					" (seq="+read_sequence.length()+", qual="+quality.length()+"). Processing read leniently.");
			}
			long estimatedReadBytes=Math.min(PROCESSING_BUDGET_BYTES,
				(2L<<20)+40L*(long)read_sequence.length());
			while (tasks_in_flight>0 &&
				(tasks_in_flight>=MAX_TASKS_IN_FLIGHT ||
				 estimatedBytesInFlight+estimatedReadBytes>PROCESSING_BUDGET_BYTES ||
				 heapUsageFraction()>0.88d))
			{
				if (tasks_in_flight>=MAX_TASKS_IN_FLIGHT) taskWindowWaits++;
				else if (estimatedBytesInFlight+estimatedReadBytes>PROCESSING_BUDGET_BYTES) processingBudgetWaits++;
				else heapPressureWaits++;
				CuttlefishResult ready=completionService.take().get();
				tasks_in_flight--;
				estimatedBytesInFlight-=ready.estimated_memory_bytes;
				reads_processed+=consumeResult(ready,panGenomeSummaries,bw,desired_confidence_or_minFreq);
			}
			completionService.submit(
				new ProcessReadTask(
					read_header,
					read_sequence,
					k,
					step,
					pp,
					nc,
					ncr,
					gct,
					bfh,
					dbdir,
					estimatedReadBytes
				)
			);
			tasks_in_flight++;
			reads_submitted++;
			estimatedBytesInFlight+=estimatedReadBytes;
			peakEstimatedBytesInFlight=Math.max(peakEstimatedBytesInFlight,estimatedBytesInFlight);
			peakTasksInFlight=Math.max(peakTasksInFlight,tasks_in_flight);
			maximumReadLength=Math.max(maximumReadLength,read_sequence.length());
			totalReadLength+=read_sequence.length();	
			if (!nextHeaderAlreadyRead)
			{
				line = readFastqLine(fastqReader, pendingFastqLines);
			}
			if (Math.random()<1/6666d)
			{
				elapsedTime = System.currentTimeMillis() - startTime;
				System.out.println("\t"+reads_processed+" reads processed in "+elapsedTime/1000+" s ("+usedram/1048576+" MB RAM used)");
			}
		}
		while(tasks_in_flight>0)
		{
			CuttlefishResult result=completionService.take().get();
			tasks_in_flight--;
			estimatedBytesInFlight-=result.estimated_memory_bytes;
			reads_processed+=consumeResult(result,panGenomeSummaries,bw,desired_confidence_or_minFreq);
		}

		elapsedTime = System.currentTimeMillis() - startTime;
		usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		System.out.println("\tAll reads ("+reads_processed+") processed in "+elapsedTime/1000+" s ("+usedram/1048576+" MB RAM used)");
		System.out.printf(Locale.ROOT,"\tAdaptive scheduling: average read length %.1f; maximum %,d; peak in-flight %,d reads / %.1f MiB%n",
			reads_processed==0?0.0:(double)totalReadLength/reads_processed,maximumReadLength,
			peakTasksInFlight,peakEstimatedBytesInFlight/1048576.0);
		System.out.printf(Locale.ROOT,"\tBackpressure events: task window %,d; processing budget %,d; heap pressure %,d%n",
			taskWindowWaits,processingBudgetWaits,heapPressureWaits);
		System.out.print("Writing pan-genome summaries on to file..");
		for(Map.Entry<Integer,CuttlefishSpecies> entry : panGenomeSummaries.entrySet())
		{
			Integer genoId = entry.getKey();
			CuttlefishSpecies geno = entry.getValue();
			//Float nkm = Float.NaN;
			//try { nkm = (float)geno.curKmerMinimiz.cardinality(); } catch (Exception e) {e.printStackTrace();}
			//if (nkm.isNaN() || nkm>0) {pw.write(genoId+","+geno.speciesName+","+Math.min(100,100f*nkm/geno.totKmerMinimiz)+"%,"+(float)geno.readDepth/nkm); pw.newLine();}
			double nkm = geno.curKmerMinimiz==null?0.0:geno.curKmerMinimiz.estimate();
			if (Double.isFinite(nkm) && nkm>0.0) {pw.write(genoId+","+csv(geno.speciesName)+","+Math.min(100.0,100.0*nkm/(double)geno.totKmerMinimiz)+"%,"+(double)geno.readDepth/nkm); pw.newLine();}
			if (Math.random()<0.05) {pw.flush();}
		}
		pw.flush();
		System.out.println(" done.");
		
		elapsedTime = System.currentTimeMillis() - timeZero;
		System.out.println("Total time employed  = "+elapsedTime/1000+" s");
		}
		catch(Exception | Error failure)
		{
			primaryFailure = failure;
			throw failure;
		}
		finally
		{
			Exception cleanupFailure = closeClassificationResources(
				fastqReader,bw,pw,bfh,executor);
			if(cleanupFailure != null)
			{
				if(primaryFailure != null) primaryFailure.addSuppressed(cleanupFailure);
				else throw cleanupFailure;
			}
		}
	}

	private static Exception closeClassificationResources(
		BufferedReader fastqReader,BufferedWriter bw,BufferedWriter pw,
		HashMap<String,ShardedBloomFilter> filters,ExecutorService executor)
	{
		Exception first = null;
		if(executor != null)
		{
			executor.shutdown();
			try
			{
				if(!executor.awaitTermination(30L,TimeUnit.SECONDS))
				{
					executor.shutdownNow();
					executor.awaitTermination(30L,TimeUnit.SECONDS);
				}
			}
			catch(InterruptedException interrupted)
			{
				executor.shutdownNow();
				Thread.currentThread().interrupt();
				first=appendCleanupFailure(first,interrupted);
			}
		}
		try { if(fastqReader != null) fastqReader.close(); }
		catch(Exception failure) { first=appendCleanupFailure(first,failure); }
		try { if(bw != null) bw.close(); }
		catch(Exception failure) { first=appendCleanupFailure(first,failure); }
		try { if(pw != null) pw.close(); }
		catch(Exception failure) { first=appendCleanupFailure(first,failure); }
		if(filters != null)
		{
			for(ShardedBloomFilter filter:filters.values())
			{
				try { if(filter != null) filter.close(); }
				catch(Exception failure) { first=appendCleanupFailure(first,failure); }
			}
		}
		return first;
	}

	private static Exception appendCleanupFailure(Exception first,Exception next)
	{
		if(first == null) return next;
		first.addSuppressed(next);
		return first;
	}
	
	private static double heapUsageFraction()
	{
		Runtime r=Runtime.getRuntime();
		return (double)(r.totalMemory()-r.freeMemory())/(double)r.maxMemory();
	}

	private static int consumeResult(CuttlefishResult cr,
		HashMap<Integer,CuttlefishSpecies> summaries,BufferedWriter bw,float threshold) throws Exception
	{
		CuttlefishSpecies cs=summaries.get(cr.species_index);
		String name=cs==null?null:cs.speciesName;
		float observed=cr.prob_or_score/(1+cr.prob_or_score);
		if (Float.isNaN(observed)) observed=0;
		if (threshold>=1) observed=cr.num_read_geno_or_cluster_hits;
		if (name==null || cr.num_read_geno_or_cluster_hits<=0 || observed<threshold) name="?";
		float match=cr.num_read_hashes<=0?0f:cr.num_read_geno_or_cluster_hits/cr.num_read_hashes;
		float confidence=cr.prob_or_score/(1+cr.prob_or_score);
		if (Float.isNaN(confidence)) confidence=0;
		bw.write(csv(cr.read_header)+","+csv(name)+","+cr.num_read_geno_or_cluster_hits+"|"+cr.num_read_hashes+","+match+","+confidence);
		bw.newLine();
		if (!name.equals("?") && cs!=null)
		{
			cs.readDepth+=cr.species_kmer_hits.size();
			if (cs.curKmerMinimiz==null) cs.curKmerMinimiz=CounterFactory.create();
			LongIterator it=cr.species_kmer_hits.iterator();
			while(it.hasNext())
			{
				long rawMinimizer=it.nextLong();
				long mixedHash=DistinctHashMixer.mix(rawMinimizer);
				cs.curKmerMinimiz.addHash(mixedHash);
			}
		}
		return 1;
	}

	public static long openFilters(HashMap<String,ShardedBloomFilter> bfh,GenomesClustersTree gct,String dbdir,long remainingBudgetBytes,AtomicInteger filtersLoadedInRAM,AtomicInteger filtersLoadedMMAP) throws Exception
	{
		if(gct==null)return remainingBudgetBytes;
		for(String bloomId:gct.clustersBloomIds)
		{
			String base=dbdir+"/"+bloomId;
			long filterBytes=ShardedBloomFilter.getDiskSizeBytesForBase(base);
			boolean loadIntoRam=filterBytes<=remainingBudgetBytes;
			ShardedBloomFilter bf=new ShardedBloomFilter(base,loadIntoRam?filterBytes:0L);
			bfh.put(bloomId,bf);
			if(loadIntoRam){remainingBudgetBytes-=filterBytes;filtersLoadedInRAM.incrementAndGet();}
			else filtersLoadedMMAP.incrementAndGet();
		}
		for(GenomesClustersTree child:gct.childrenList)
			remainingBudgetBytes=openFilters(bfh,child,dbdir,remainingBudgetBytes,filtersLoadedInRAM,filtersLoadedMMAP);
		return remainingBudgetBytes;
	}

	public static float[] matchScores(int queryElements, ArrayList<Integer> clusterElements, IntArrayList queryHits) throws Exception
	{
		if(queryElements<0 || clusterElements==null || queryHits==null || clusterElements.isEmpty() || clusterElements.size()!=queryHits.size()) throw new IllegalArgumentException("Invalid scoring inputs");
		long totalClusterElements = 0L;
		for (int n : clusterElements) {if(n<=0)throw new IllegalArgumentException("Invalid zero or negative cluster size"); totalClusterElements=Math.addExact(totalClusterElements,(long)n);}
		if(queryElements==0 || totalClusterElements==0L)return new float[clusterElements.size()];
		float[] scores = new float[clusterElements.size()];
		for (int i = 0; i < clusterElements.size(); i++)
		{
			int clusterSize = clusterElements.get(i);
			int hits = queryHits.getInt(i);
			if(hits<0)throw new IllegalArgumentException("Negative hit count");
			double expected = (double) queryElements * clusterSize / totalClusterElements;
			double score;
			
			//unsigned deviance
			//if (hits <= expected) {score = 0.0;} else {score = 2.0 * (hits * Math.log((double) hits / expected) - (hits - expected));} 
			
			//signed deviance
			//if (hits == 0) {score = -2.0 * expected;} else {double deviance = 2.0 * (hits * Math.log((double) hits / expected) - (hits - expected)); score = Math.signum(hits - expected) * Math.sqrt(deviance);}
			
			//hybrid sqrt(M)R
			double enrichment = (double) hits / expected; score = Math.sqrt(hits) * enrichment;
			
			scores[i] = (float) score;
		}
		return scores;
	}
	
	public static float cumulProba(double x, double m, double s)
	{
		NormalDistribution normDist = new NormalDistribution(m,s);
		double pValu = normDist.cumulativeProbability(x-1);
		//PoissonDistribution poisDist = new PoissonDistribution(m);
		//double pValu = poisDist.cumulativeProbability((int)(x-1));
		return (float)pValu;
	}


	private static String readFastqLine(
		BufferedReader br,
		ArrayDeque<String> pendingLines
	) throws IOException
	{
		return pendingLines.isEmpty() ? br.readLine() : pendingLines.removeFirst();
	}

	private static boolean isPlausibleNextFastqHeader(
		BufferedReader br,
		ArrayDeque<String> pendingLines
	) throws IOException
	{
		ArrayList<String> lookahead = new ArrayList<String>();
		boolean sawSequence = false;
		boolean plausible = false;
		final int MAX_LOOKAHEAD_LINES = 100000;

		for (int i=0; i<MAX_LOOKAHEAD_LINES; i++)
		{
			String candidate = readFastqLine(br, pendingLines);
			if (candidate == null)
			{
				break;
			}
			lookahead.add(candidate);
			if (candidate.startsWith("+"))
			{
				plausible = sawSequence;
				break;
			}
			if (candidate.startsWith("@") || !isPlausibleSequenceLine(candidate))
			{
				break;
			}
			sawSequence = true;
		}

		for (int i=lookahead.size()-1; i>=0; i--)
		{
			pendingLines.addFirst(lookahead.get(i));
		}
		return plausible;
	}

	private static boolean isPlausibleSequenceLine(String line)
	{
		if (line.length() == 0) return false;
		for (int i=0; i<line.length(); i++)
		{
			char c = Character.toUpperCase(line.charAt(i));
			if ("ACGTURYSWKMBDHVN.-*".indexOf(c) < 0)
			{
				return false;
			}
		}
		return true;
	}

	private static final String EXPECTED_DB_FORMAT="2";
	private static final String EXPECTED_MINIMIZER="canonical-2bit-exact-v1";
	private static final String EXPECTED_MINIMIZER_WINDOW="k-step-adjacent-mmers-v1";
	private static final String EXPECTED_BLOOM_FORMAT="legacy-header-prime-bits-v1";
	private static final String EXPECTED_BLOOM_HASH="murmurhash2-32-seeds-0-1-big-endian-long";
	private static final String EXPECTED_BLOOM_PROBE="double-hash-prime-modulus-zero-step-to-one-v1";
	private static final String EXPECTED_BLOOM_SHARD_SELECTION="seed0-low-bits-power-of-two-v1";
	private static Properties loadAndValidateMetadata(String dbdir) throws IOException
	{
		Path path=Paths.get(dbdir,"metadata.txt");
		if(!Files.isRegularFile(path))throw new IOException("Missing database metadata: "+path);
		Properties p=new Properties();
		try(BufferedReader reader=Files.newBufferedReader(path,StandardCharsets.UTF_8)){p.load(reader);}
		if(!"true".equals(p.getProperty("buildComplete")))throw new IOException("Database build is incomplete: "+dbdir);
		requireMetadataValue(p,"cuttlefishDbFormat",EXPECTED_DB_FORMAT);
		requireMetadataValue(p,"minimizerAlgorithm",EXPECTED_MINIMIZER);
		requireMetadataValue(p,"minimizerWindowDefinition",EXPECTED_MINIMIZER_WINDOW);
		requireMetadataValue(p,"bloomFormat",EXPECTED_BLOOM_FORMAT);
		requireMetadataValue(p,"bloomHashAlgorithm",EXPECTED_BLOOM_HASH);
		requireMetadataValue(p,"bloomProbePolicy",EXPECTED_BLOOM_PROBE);
		requireMetadataValue(p,"bloomShardSelection",EXPECTED_BLOOM_SHARD_SELECTION);
		for(String required:new String[]{"taxonomicClusterTree.ser","info.txt"})
			if(!Files.isRegularFile(Paths.get(dbdir,required)))throw new IOException("Missing database file: "+required);
		return p;
	}

	private static void requireMetadataValue(Properties p,String key,String expected) throws IOException
	{
		String actual=p.getProperty(key);
		if(!expected.equals(actual))
			throw new IOException("Incompatible database metadata "+key+": expected "+expected+", found "+actual);
	}
	private static int parsePositiveMetadataInt(Properties p,String key) throws IOException
	{
		String value=p.getProperty(key);
		if(value==null)throw new IOException("Missing database metadata value: "+key);
		try
		{
			int parsed=Integer.parseInt(value);
			if(parsed<=0)throw new IOException("Invalid database metadata "+key+": must be positive");
			return parsed;
		}
		catch(NumberFormatException invalid)
		{
			throw new IOException("Invalid integer database metadata "+key+": "+value,invalid);
		}
	}

	private static String csv(String value)
	{
		if(value==null)return "";
		if(value.indexOf(',')<0&&value.indexOf('"')<0&&value.indexOf('\r')<0&&value.indexOf('\n')<0)return value;
		return "\""+value.replace("\"","\"\"")+"\"";
	}
	private static String requireValue(String[] args, int index, String option)
	{
		if (index >= args.length)
		{
			throw new IllegalArgumentException(
				"Missing value after " + option
			);
		}
		return args[index];
	}
	
}


final class DistinctHashMixer
{
	private static final long SEED=0x9e3779b97f4a7c15L;
	private DistinctHashMixer(){}
	static long mix(long minimizer)
	{
		long x=minimizer^SEED;
		x^=x>>>33;
		x*=0xff51afd7ed558ccdL;
		x^=x>>>33;
		x*=0xc4ceb9fe1a85ec53L;
		x^=x>>>33;
		return x;
	}
}

interface DistinctCounter
{
	void addHash(long hash);
	double estimate();
}
final class UltraDistinctCounter implements DistinctCounter
{
	private final UltraLogLog sketch;
	UltraDistinctCounter(int precision){sketch=UltraLogLog.create(precision);}
	public void addHash(long hash){sketch.add(hash);}
	public double estimate(){return sketch.getDistinctCountEstimate();}
}
final class CompactDistinctCounter implements DistinctCounter
{
	private final HllSketch sketch;
	CompactDistinctCounter(int precision){sketch=new HllSketch(precision,TgtHllType.HLL_4);}
	public void addHash(long hash){sketch.update(hash);}
	public double estimate(){return sketch.getEstimate();}
}
final class CounterFactory
{
	private static boolean ultra=true;
	private static int precision=12;
	static void configure(boolean useUltra,int p){ultra=useUltra; precision=p;}
	static DistinctCounter create(){return ultra?new UltraDistinctCounter(precision):new CompactDistinctCounter(precision);}
}

class CuttlefishSpecies implements Serializable
{
	public String speciesName;
	public long totKmerMinimiz;
	public DistinctCounter curKmerMinimiz;
	public long readDepth;
	public CuttlefishSpecies()
	{
		speciesName = "";
		totKmerMinimiz = 0;
		curKmerMinimiz = null;
		readDepth = 0;
	}
}

class CuttlefishResult implements Serializable
{
	public String read_header;
	public int read_sequence_length;
	public LongOpenHashSet read_hashes;
	public long[] read_hash_array;
	public int num_read_hashes;
	public long estimated_memory_bytes;
	public int species_index;
	public float num_read_geno_or_cluster_hits;
	public float prob_or_score;
	public LongArrayList species_kmer_hits;
	
	public CuttlefishResult()
	{
		read_header=null;
		read_sequence_length=-1;
		read_hashes=new LongOpenHashSet();
		read_hash_array=null;
		num_read_hashes=0;
		estimated_memory_bytes=0L;
		species_index=-1;
		num_read_geno_or_cluster_hits=-1f;
		prob_or_score=-3.4E38f;
		species_kmer_hits=new LongArrayList();
	}
	
	public boolean classifySequence(GenomesClustersTree gct, HashMap<String,ShardedBloomFilter> bfh, String dbdir) throws Exception
	{
		IntArrayList numReadHits = new IntArrayList();
		ArrayList<LongArrayList> kmerReadHits = new ArrayList<LongArrayList>();
		for (int i=0; i<gct.clustersIds.size(); i++)
		{
			int localScore = 0;
			LongArrayList localKmers = new LongArrayList();
			boolean bloomFilterInCache = true;
			ShardedBloomFilter bf = bfh.get(gct.clustersBloomIds.get(i));
			if (bf==null)
			{
				bf = new ShardedBloomFilter(dbdir+"/"+gct.clustersBloomIds.get(i),0);
				bloomFilterInCache = false;
			}
			try
			{
				long[] hashes=this.read_hash_array;
				for (int h=0; h<hashes.length; h++)
				{
					long kmerHash=hashes[h];
					if (bf.contains(kmerHash)) {localScore++; localKmers.add(kmerHash);}
				}
			}
			finally
			{
				if (!bloomFilterInCache) bf.close();
			}
			numReadHits.add(localScore);
			kmerReadHits.add(localKmers);
		}
		float [] zvals = CUTTLEFISH.matchScores(this.read_hashes.size(),gct.clusterElements,numReadHits);
		int bestIndex = 0;
		for (int i=1; i<zvals.length; i++) {if (zvals[i]>zvals[bestIndex]) {bestIndex=i;}}
		int bestScore = numReadHits.getInt(bestIndex);
		this.species_index = gct.clustersIds.get(bestIndex).get(0);
		this.num_read_geno_or_cluster_hits = bestScore;
		this.prob_or_score = zvals[bestIndex]; 
		if (bestScore<=0) 
		{
			return false;
		}
		if (bestScore>0 && gct.clustersIds.get(bestIndex).size()==1)
		{
			this.species_kmer_hits = kmerReadHits.get(bestIndex);
			return true;
		}
		return classifySequence(gct.childrenList.get(bestIndex),bfh,dbdir);
	}
}

class ProcessReadTask implements Callable<CuttlefishResult>
	{
		private String rh;
		private String rs;
		private int k;
		private int step;
		private long[][] pp;
		private int[] nc;
		private int[] ncr;
		private GenomesClustersTree gct;
		private HashMap<String,ShardedBloomFilter> bfh;
		private String dbdir;
		private long estimatedMemoryBytes;
		public ProcessReadTask()
		{		
			rh=null;
			rs=null;
			k=0;
			step=0;
			pp=null;
			nc=null;
			ncr=null;
			gct=null;
			bfh=null;
			dbdir=null;
		}
		public ProcessReadTask(String readhdr, String readseq, int kappa, int stepminimiz, long[][] precpows, int[] nch, int[] ncrh, GenomesClustersTree genclutre, HashMap<String,ShardedBloomFilter> blofilhas, String dbdire, long estimatedBytes)
		{
			rh=readhdr;
			rs=readseq;
			k=kappa;
			step=stepminimiz;
			pp=precpows;
			nc=nch;
			ncr=ncrh;
			gct=genclutre;
			bfh=blofilhas;
			dbdir=dbdire;
			estimatedMemoryBytes=estimatedBytes;
		}
		public CuttlefishResult call() throws Exception 
		{
			CuttlefishResult result = new CuttlefishResult();
			result.read_header = rh;
			result.read_sequence_length = rs.length();
			result.read_hashes = new LongOpenHashSet();
			BuildCUTTLEFISH.hashSequenceInto(rs,k,step,pp,nc,ncr,result.read_hashes);
			result.num_read_hashes=result.read_hashes.size();
			result.read_hash_array=result.read_hashes.toLongArray();
			boolean classified = result.classifySequence(gct,bfh,dbdir);
			result.read_hashes=null;
			result.read_hash_array=null;
			result.estimated_memory_bytes=estimatedMemoryBytes;
			return result;
		}
	}


