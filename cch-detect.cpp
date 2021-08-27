// Brief background
// ----------------
// 
// Given a biallelic locus with alleles A and B, and two samples, one homozygous AA
// and one homozygous BB, we say this locus exhibits conflicting homozygozity (CH).
// Short of a mutation, this means that genetic material spanning that locus can not
// be identical by descent (IBD).
//
// Conversely, a long run of consecutive loci with no conflicting homozygozity (a
// CH-free region) is suggestive of IBD.
// 
// Across a larger number of samples N, if such a large (4 cM or longer) region occurs
// in all samples, it very strongly suggests IBD in that region. If such a region
// exists only in a subset of K samples, it follows that it is only IBD across the K
// samples. Even if all N samples share a phenotype, it likely has a different cause
// in the other N - K samples. This allows the identification of a region of interest
// even in the face of such confounding phenotypes.

// Algorithm
// ---------
// 
// This program finds long (> 4 cM) regions of no-CH across each combination of K
// samples out of the given N.
// 
// The program looks at all N-choose-K combinations, but delivers reasonable speed
// thanks to a number of techniques that, sadly, make it harder to read.
// 
// 1. Data representation
// 
// For a given sample, a genotype is represented as one bit in each of two bit vectors
// (arrays of 64-bit integers), named AA and BB. The bit is set in AA if the sample is
// homozygous AA, and the corresponding bit is set in BB if the sample is homozygous
// BB. The same representation works for multiple samples, where the bit in AA is set
// if any of the samples are homozygous AA, and likewise BB. The bit vectors for a set
// of samples are the bitwise-OR of the bit vectors for the individual samples --
// crucially, an associative operation (and fast, 64-way parallel on a 64-bit CPU).
// 
// Clearly, for a single sample, the same bit is never set in both the AA and BB bit 
// vectors (the same locus cannot be both homozygous AA and homozygous BB). For
// multiple samples, a bit set in both AA and BB marks conflicting homozygozity; this
// is checked with a bitwise-AND, 64 loci at a time. If there is CH, this is non-zero.
// As described below, only the position of the first or last conflict is needed,
// which is found quickly with the LZCNT and TZCNT instructions (leading zero count,
// trailing zero count) available on modern CPUs.
// 
// 2. Slices
// 
// A key observation is that, if the chromosomes are divided into 4 cM slices, a
// region long enough to be interesting (4 cM or more) will cross (or, at least,
// touch) a slice boundary. It can not be wholly contained within a slice, simply
// because it's larger than a slice. Thus, it must begin in one slice and end in a
// different slice (spanning zero or more other slices in between).
// 
// We can thus reduce the information about a slice and a given combination of samples
// to two numbers: p, the distance from beginning of the slice to the first CH, and q,
// the distance from the last CH to the end of the slice -- plus a special case, when
// there is no CH within the slice.
// 
// As we process the data one slice at a time, we need to retain only one number per
// combination, namely the distance from the last encountered CH to the end of the
// last processed slice. Equivalently, this is the length of the current CH-free run,
// up to the slice boundary. This is 'chr_q', which is the 'q' of the last slice with
// a CH, plus the length of any CH-free slices following it. If the current slice is
// CH-free, it extends this CH-free run to the next slice boundary. If it's not, its
// first CH ends the run (run length: chr_q + p) and its last CH becomes the start of
// the new current run.
// 
// Regions do not span chromosome boundaries, so there's logic to terminate the slice
// and end current regions at the end of a chromosome, then re-initialize the
// algorithm for the next chromosome in the input file.
//
// The actual logic takes into account the number of SNPs in a region in addition to
// its genetic length, and works mostly with the indices of SNPs rather than directly
// with genetic distances, but the principle is as outlined above.
//
// Memory usage does not depend on the number of SNPs in the input file. The file is
// read only once, sequentially.
// 
// 3. Recursive generation of combinations
// 
// The N-choose-K combinations of samples are computed using what we believe to be a
// novel recursive approach. Each combination is the union of two smaller
// combinations, each with K/2 elements (plus one extra element, if K is odd).
// This is not a particularly effective way of generating combinations if all one
// needs is a list, but it's very suitable for computing a property of all
// combinations if it can easily be computed from the properties of the two halves.
// 
// Putting together two subsets of samples will never make a CH-free region larger,
// only smaller. The CH-free regions of the resulting set are contained within the
// CH-free regions either subset -- they're contained within the intersection of the
// CH-free regions of the two subsets. This is obvious, because adding more samples
// can only create new conflicts, not remove them.
// 
// In the context of per-slice processing, this means that, if either half of the
// combination does not have a region long enough, their union does not. If they do,
// then any new conflicts are found by OR-ing the AA bit vectors of the two halves,
// OR-ing their BB bit vectors, and AND-ing the two results. This calculation is done
// only within the intersection of the regions of the two combination halves.
// 
// In the case of the "final" combinations (K elements, not used to make larger
// combinations -- stored in 'final_combs'), the bit vectors don't need to be stored,
// only calculated once, and only where not ruled out as above. In the case of the
// smaller combinations that make up larger ones (in the vector 'combs'), the bitmaps
// are stored (in 'bitmaps'), because each may be used multiple times.
// 
// A further optimization avoids some of the OR-ing and AND-ing by keeping the bitmaps
// of both the current and the previous slice. We defer the finding of the last CH in
// the previous slice until after we find the first CH in the current one, which can
// obviate the need if it leads to too short a region anyway. This doubles the size of
// 'bitmaps' but reduces execution time.
// 
// A combination, regardless of size, is stored as three numbers: the indices X and Y
// of its two halves, and, if needed, the index S of the extra element. The sets of
// elements that make each combination are not stored explicitly; they are
// reconstructed from this tree when needed (i.e. for printing out). As stated above,
// the algorithm needs only one other piece of data per combination, so, at 16 bytes
// per combination, we can hold many of them in memory at the same time.
//
// 4. Integer arithmetic
//
// Genetic positions are scaled to micromorgans and stored as integers. The lengths
// of all human chromosomes fit into 32-bit integers, and four decimals of precision
// for the value in cM is enough.

// Usage
// -----
//
// This is not a finished program and has no command line processing yet. Search for
// "TODO command line" comments.
//
// N, the number of samples, is deduced from the first line of the input file.
// K, the size of the combinations, is hardcoded, but 0 means use K = N, and
//    a negative value means use |K| less than N.
// min_gen_dist, default 4 cM, the minimum size of the CH-free regions to detect;
//    larger is faster!
// min_snp_count, default 100, the minimum number of SNPs in a region.
// max_slice_snp_count, default 640, multiple of 64, the maximum size of a slice.
//    Choose according to genotyping density, so that almost all regions of the chosen
//    size have fewer SNPs than this. Larger uses more memory; too small makes the
//    slicing ineffective and hurts performance. The default works well for 4 cM and
//    270,000 SNPs. (Multiply the region size threshold in cM by the number of SNPs in
//    the whole genome, divide by 1700, round up.) This should remain a compile-time
//    parameter; making it variable sacrifices some performance.
//
// The input file name is hardcoded in main(), at the bottom of this file.
// Output is written to standard output.
// Progress information is written to standard error.
//
// The input file is a tab-delimited text file, without any comments or header
// (no sample IDs). Each line must contain:
// - chromosome number, an integer (no 'chr')
// - marker name, text (not used)
// - genetic position, a real number, in cM
// - physical position, an integer (not important whether 0-based or 1-based)
// - exactly N genotypes, which can be 'AA', 'AB', or 'BB' (but not 'BA').
// The lines must be grouped by chromosome and sorted by position within a chromosome.
// Chromosome order is not important (so 1, 10, 11, ..., 19, 2, 20, 21, 22 is fine).
// If you have X and Y chromosomes, use a number (e.g. 23, 24).
//
// If roughly the same region occurs in K+1 combinations, increase K by 1.
// If it occurs in (K+1)(K+2) combinations, increase K by 2. In fact, there's usually
// no need to rerun the program -- the intersection of the printed regions is the
// correct result for the larger K, and the combination it applies to is the union of
// the printed combinations.
// TODO This hints at a method to find maximal subsets without computing all
// TODO    combinations, which may be faster and use less memory.
//
// If there is no output, decrease K. Maybe try a small number (say, 7) and
// interrupt if gigantic output.

// Limitations
// -----------
//
// Execution time increases dramatically with both N and K.
//
// All combinations must fit in memory, and the program is single-threaded -- but
// see the musings after the class definition.
//
// The program does not support filtering or subsetting the samples or the SNPs.
//
// The program does not have a feature to tolerate a small number (e.g. 1 or 2) of
// conflicts within a region. The same overall structure works with the obvious
// dynamic programming solution, at the cost of making p, q, and chr_q arrays, and
// complicating the per-slice logic to find not just the first and last CH locus.

// Compilation
// -----------
//
// On Windows, compile as 64-bit application.
//
// On Linux and other gcc targets, use
//    g++ -march=x86-64 -mlzcnt -mbmi -O3 cch-detect.cpp -o cch-detect


#include <string>
#include <vector>
#include <array>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#if defined(_WIN32)
#include <intrin.h>
#elif defined(linux)
#include <x86intrin.h>  // needs g++ -march=x86-64 -mlzcnt -mbmi (and don't forget -O3)
#else
#error("Don't know how to get the LZCNT and TZCNT intrinsics")
#endif


#define GEN_DIST_SCALE 10000.0  // use micromorgans for genetic distances internally; 1 cM = 10000 uM

class Detector
{
public:
	typedef int chr_t;
	typedef int gen_pos_t;
	typedef int phys_pos_t;

	typedef unsigned long long bitmap_word_t;
	static constexpr int snp_word_bits = 64;  // size, in bits, of bitmap_word_t, and a power of 2

public:
	gen_pos_t min_gen_dist = static_cast<gen_pos_t>(4.0 * GEN_DIST_SCALE);  // 4 cM // TODO command line
	int min_snp_count = 100;  // TODO command line
	static constexpr int max_slice_snp_count = 640;  // best if a multiple of snp_word_bits
	static constexpr int max_slice_word_count = (max_slice_snp_count + snp_word_bits - 1) / snp_word_bits;

	struct comb_t
	{
		// X, S and Y are indices into the vector of combinations.
		// If S == -1, this combination is the union of the combinations X and Y;
		// if S >= 0, that of X, S and Y. "Combination" S is a single element, so S < N.
		// For single-element "combinations", X = Y = -1.
		int X;
		int S;
		int Y;
	};

	struct final_comb_t
	{
		int X;
		int S;
		int Y;
		int chr_q;  // number of consecutive no-conflict SNPs just before the current slice
		            // or -1 if not computed, but falls somewhere within the previous slice
	};

	struct bitmap_comb_t
	{
		// One bit per SNP, LSB first, set if any member of the combination is...
		bitmap_word_t AA[max_slice_word_count];  // ... homozygous AA
		bitmap_word_t BB[max_slice_word_count];  // ... homozygous BB
		int p;  // number of consecutive no-conflict SNPs at the beginning of the slice
		int q;  // and at the end of the slice
	};

	struct snp_t
	{
		std::string name;
		gen_pos_t gen_pos;
		phys_pos_t phys_pos;
	};

	int N = 0;
	int K = 9;  // TODO command line
	std::vector<comb_t> combs;
	std::vector<final_comb_t> final_combs;
	std::vector<std::array<bitmap_comb_t, 2>> bitmaps;
private:
	std::vector<int> next;
	std::map<std::tuple<int, int, int>, int> mem;

public:
	chr_t chr = 0;
	int crt_slice = 0;
	gen_pos_t gen_pos = 0;
	gen_pos_t prev_slice_gen_pos = 0;
	int slice_snp_count = 0;
	int prev_slice_snp_count = 0;
	int slice_word_count = 0;
	int prev_slice_word_count = 0;
	int slice_padding = 0;
	int prev_slice_padding = 0;
	bitmap_word_t bit_mask = 1;
	int slice_start_snp = 0;
	int min_range_snp = 0;
	std::vector<int> min_range_q;
	std::vector<snp_t> chr_snps;

public:
	void process_input(std::istream& in);
	void compute_combinations();
	int comb(int x, int y, int k, bool is_final);
	void add_comb(int x, int s, int y, int& head, int& tail, bool is_final);
	void extract_comb(int comb_index, std::vector<int>& C);
	void extract_final_comb(const final_comb_t& comb, std::vector<int>& C);
	void begin_slice();
	void next_snp();
	void end_slice(bool is_last);
	void begin_chromosome();
	void end_chromosome();

private:
	void compute_homozygosity_bitmaps(int comb_index);
	int extend_conflict_free_region(final_comb_t comb, bool chromosome_ends);
};


// Multithreading notes
// --------------------
//
// If dividing work by chromosome:
// The combinations (combs, and final_comb_t except chr_q) can be shared by all
// threads. Each thread will need its own copy of bitmaps, chr_snps, min_range_q etc,
// and final_comb_t::chr_q would need to go in a separate per-thread structure.
// This would need extra memory.
//
// If dividing the work by combination:
// Once combs are computed, a pool of threads can calculate final_combs in any
// order, because they're independent. They'd need to rendezvous to continue to
// the next slice. Because combs depend on each other, multithreading that
// computation would require some synchronization, or some reordering (e.g. by
// combination size, because those of the same size are independent, so they can
// be handled by a pool of threads). If ordered as above, just provide an atomic
// index of 'next comb available for processing'; an available thread takes on
// that comb and increments the index. Also add a 'completed' flag to each
// comb. If the prerequisites are not both 'completed', a thread can wait and
// retry -- or maybe add itself to a list of threads to be woken up when specific
// combs get completed. This might work well enough even with unordered combs.
// Or order the combs, and make sure all threads are done with combinations of
// one size before moving on to the next size; no 'completed' flag is needed then.


// More combinations than fit in memory
// ------------------------------------
//
// The program can be modified to calculate only some combinations, not all of
// them, for instance by stopping generating final_combs after a certain number
// is reached. It could then be run again with the next batch of combinations,
// because the state of the combination generator can be saved. At a minimum,
// it would need the x and y from compute_combinations(), and a count of how
// many have already been generated from that series, so they can be skipped.


// SNP data processing
// ----------------------------------------------------------------------------

void Detector::process_input(std::istream& in)
{
	chr = 0;
	std::string line;

	while (std::getline(in, line))
	{
		// TODO read header and extract individual IDs
		// TODO skip comment lines
		// TODO validate the data, at least a little

		auto sep = line.find('\t');
		chr_t new_chr = std::stoi(line.substr(0, sep));

		auto pos = sep + 1;
		sep = line.find('\t', pos);
		std::string marker_name = line.substr(pos, sep - pos);

		pos = sep + 1;
		sep = line.find('\t', pos);
		double unscaled_gen_pos = std::stod(line.substr(pos, sep - pos));
		gen_pos_t new_gen_pos = std::lround(unscaled_gen_pos * GEN_DIST_SCALE);

		pos = sep + 1;
		sep = line.find('\t', pos);
		phys_pos_t phys_pos = std::stoi(line.substr(pos, sep - pos));

		int new_N = (static_cast<int>(line.length() - sep)) / 3;

		if (N == 0)
		{
			N = new_N;
			if (K <= 0)
				K = N + K;  // 0 means same as N, negative means |K| less than N
			K = std::min(N, std::max(4, K));
			compute_combinations();
		}

		if (new_chr != chr)
		{
			end_chromosome();

			chr = new_chr;
			begin_chromosome();
		}

		if (slice_snp_count == 0)
			slice_start_snp = static_cast<int>(chr_snps.size());

		gen_pos = new_gen_pos;

		chr_snps.push_back({ marker_name, gen_pos, phys_pos });

		while (gen_pos - chr_snps[min_range_snp].gen_pos >= min_gen_dist
			&& chr_snps.size() - min_range_snp >= min_snp_count)
			++min_range_snp;

		min_range_q.push_back(slice_start_snp - min_range_snp + 1);

		for (int i = 0; i < N; ++i, sep += 3)
		{
			// The rest of the line consists of three-character groups, one per
			// individual. There are only three possibilities: "\tAA", "\tAB", "\tBB"
			// ("\tBA" is not allowed).

			// If homozygous, set the respective bit in either the AA bitmap or the BB bitmap.
			// If heterozygous, leave both bits zero.
			if (line[sep + 2] == 'A')  // '?A' can only be 'AA'
				bitmaps[i][crt_slice].AA[slice_word_count] |= bit_mask;
			else if (line[sep + 1] == 'B')  // 'B?' can only be 'BB'
				bitmaps[i][crt_slice].BB[slice_word_count] |= bit_mask;
		}

		next_snp();
	}

	end_chromosome();

	std::cerr << "Done.                               " << std::endl;  // spaces wipe progress message
}

void Detector::begin_chromosome()
{
	prev_slice_gen_pos = 0;
	prev_slice_snp_count = 0;
	prev_slice_word_count = 0;
	min_range_snp = 0;
	chr_snps.clear();
	for (auto& comb : final_combs)
		comb.chr_q = 0;  // this also ensures the code does not try to look at the previous slice
	crt_slice = 0;

	begin_slice();
}

void Detector::begin_slice()
{
	slice_snp_count = 0;
	slice_word_count = 0;
	bit_mask = 1;
	for (int i = 0; i < N; ++i)
		bitmaps[i] = { 0 };

	min_range_q.clear();
}

void Detector::next_snp()
{
	++slice_snp_count;
	
	// Is the slice full, or is it long enough that it can contain a region
	// of the given minimum genetic distance and minimum number of markers?
	//
	// The key to the algorithm is that any such region must cross, or at least
	// touch, a slice boundary. If a slice cannot contain an embedded region of
	// interest, we can reduce the slice to just two numbers (per combination):
	// - the number of consecutive no-conflict markers at the beginning of the
	//   slice. These extend and terminate the trailing no-conflict region
	//   from the previous slice;
	// - the number of consecutive no-conflict markers at the end of the slice.
	//   These begin a new no-conflict region that may be extended by the next
	//   slice;
	// plus a special case:
	// - no conflict throughout the slice. This extends the (possibly zero-
	//   length) no-conflict region from the previous slice, and can be further
	//   extended by the next slice.
	//
	if (slice_snp_count >= max_slice_snp_count ||
		gen_pos - prev_slice_gen_pos >= min_gen_dist && slice_snp_count >= min_snp_count)
	{
		end_slice(false);
		begin_slice();
	}
	else
	{
		bit_mask <<= 1;
		if (bit_mask == 0)
		{
			++slice_word_count;
			bit_mask = 1;
		}
	}
}

void Detector::end_slice(bool chromosome_ends)
{
	if (bit_mask != 1)
		++slice_word_count;

	// Unused bits in the last word: difference between slice_snp_count and the next
	// multiple of snp_word_bits (the latter is a power of 2)
	slice_padding = (-slice_snp_count) & (snp_word_bits - 1);

	std::cerr << "Chromosome " << std::setw(2) << std::setprecision(2) << chr << ": "
		<< std::fixed << std::setw(6) << std::setprecision(2) << (prev_slice_gen_pos / GEN_DIST_SCALE)
		<< " to " << std::setw(6) << std::setprecision(2) << (gen_pos / GEN_DIST_SCALE) << " cM  \r";

	// Add an extra element if at the end of chromosome, because the last slice is
	// handled as if there was a conflict one past the end of the last slice.
	// The extra element avoids a special case.
	if (chromosome_ends)
		min_range_q.push_back(slice_start_snp - min_range_snp + 1);

	// Single individuals: there are no conflicts
	for (int i = 0; i < N; ++i)
	{
		auto& bitmap = bitmaps[i][crt_slice];
		bitmap.p = bitmap.q = slice_snp_count;
	}

	// The small combinations used to make larger combinations
	for (int i = N; i != combs.size(); ++i)
		compute_homozygosity_bitmaps(i);

	// The final (K-element) combinations we're interested in
	for (auto& comb: final_combs)
		comb.chr_q = extend_conflict_free_region(comb, chromosome_ends);

	// Move on to the next slice
	crt_slice = !crt_slice;  // flip bitmaps: what was current slice is now previous, and vice-versa
	prev_slice_snp_count = slice_snp_count;
	prev_slice_word_count = slice_word_count;
	prev_slice_padding = slice_padding;
	prev_slice_gen_pos = gen_pos;
}

void Detector::end_chromosome()
{
	if (chr == 0)
		return;

	end_slice(true);
}


// Per-combination, per-slice calculations
// ----------------------------------------------------------------------------

void Detector::compute_homozygosity_bitmaps(int comb_index)
{
	auto& comb = combs[comb_index];
	auto& bitmap = bitmaps[comb_index][crt_slice];
	auto& X = bitmaps[comb.X][crt_slice];
	auto& Y = bitmaps[comb.Y][crt_slice];
	auto& S = bitmaps[comb.S >= 0 ? comb.S : comb.Y][crt_slice];  // TODO write the code twice, for speed gain

	bitmap.p = bitmap.q = slice_snp_count;

	// Find the first conflict and compute the bitmaps from the beginning to that point
	for (int j = 0; j < slice_word_count; ++j)
	{
		bitmap_word_t AA = X.AA[j] | Y.AA[j] | S.AA[j];
		bitmap_word_t BB = X.BB[j] | Y.BB[j] | S.BB[j];
		bitmap.AA[j] = AA;
		bitmap.BB[j] = BB;
		bitmap_word_t conflict_mask = AA & BB;
		if (conflict_mask != 0)
		{
			// Number of consecutive no-conflict SNPs at the beginning of the slice:
			// j whole words, plus the number of trailing zeros in the partial word
			// (per-SNP bits are stored LSB first, so count *trailing* zeros)
			bitmap.p = j * snp_word_bits + static_cast<int>(_tzcnt_u64(conflict_mask));

			// Find the last conflict and compute bitmaps from there to the end (work backwards)
			for (j = slice_word_count - 1; j >= 0; --j)
			{
				AA = X.AA[j] | Y.AA[j] | S.AA[j];
				BB = X.BB[j] | Y.BB[j] | S.BB[j];
				bitmap.AA[j] = AA;
				bitmap.BB[j] = BB;
				conflict_mask = AA & BB;
				if (conflict_mask != 0)
				{
					// Number of consecutive no-conflict SNPs at the end of the slice:
					// a number of whole words, plus leading zeros in the current, partial word
					// (per-SNP bits are stored LSB first, so count *leading* zeros)
					bitmap.q = (slice_word_count - 1 - j) * snp_word_bits
						+ static_cast<int>(_lzcnt_u64(conflict_mask)) - slice_padding;

					return;  // found at least one conflict
				}
			}
		}
	}

	// no conflict; bitmap.p and bitmap.q remain slice_word_count
}

int Detector::extend_conflict_free_region(final_comb_t comb, bool chromosome_ends)
{
	auto& X = bitmaps[comb.X][crt_slice];
	auto& Y = bitmaps[comb.Y][crt_slice];

	int p = std::min(X.p, Y.p);  // Upper bound for p
		// No need to use S: there's no conflict within a single individual.
		// Note p is slice_snp_count only if both halves are conflict-free, otherwise it's less

	int prev_q = comb.chr_q;
	if (prev_q == -1)  // Actual prev_q has not been calculated
	{
		auto& prev_X = bitmaps[comb.X][!crt_slice];
		auto& prev_Y = bitmaps[comb.Y][!crt_slice];
		prev_q = std::min(prev_X.q, prev_Y.q);  // Upper bound for prev_q
	}

	// If this slice definitely has a conflict (or is the last slice of a chromosome),
	// and the upper bounds for prev_q and p do not permit a region long enough,
	// there's no change of finding a long-enough region here. We're done.
	// (Note that prev_q may be the actual value or just an upper bound; this is
	// addressed further down.)
	if ((p < slice_snp_count || chromosome_ends) && prev_q < min_range_q[p])
		return -1;  // not calculated, starts in this slice

	// TODO write the code twice, for speed gain
	auto& S = bitmaps[comb.S >= 0 ? comb.S : comb.Y][crt_slice];

	// Find actual conflict, if any
	for (int j = 0; j < slice_word_count; ++j)
	{
		bitmap_word_t AA = X.AA[j] | Y.AA[j] | S.AA[j];
		bitmap_word_t BB = X.BB[j] | Y.BB[j] | S.BB[j];
		bitmap_word_t conflict_mask = AA & BB;
		if (conflict_mask != 0)
		{
			// Number of consecutive no-conflict SNPs at the beginning of the slice
			// (see similar code above)
			p = j * snp_word_bits + static_cast<int>(_tzcnt_u64(conflict_mask));
			break;
		}
	}

	// Same check, but with actual value for p: if there's a conflict in this slice and p and
	// the upper bound for prev_q do not make for a long enough region, there's nothing here.
	if ((p < slice_snp_count || chromosome_ends) && prev_q < min_range_q[p])
		return -1;  // not calculated, starts in this slice

	// If we used an upper bound for prev_q because comb.chr_q had not been calculated for
	// the previous slice, calculate it now. We waited until the last possible time
	if (comb.chr_q == -1)
	{
		auto& prev_X = bitmaps[comb.X][!crt_slice];
		auto& prev_Y = bitmaps[comb.Y][!crt_slice];
		auto& prev_S = bitmaps[comb.S >= 0 ? comb.S : comb.Y][!crt_slice];

		// Find the last conflict in the previous slice.
		// We know a previous slice exists and has a conflict, because:
		// - this is not the first slice of a chromosome (or comb.chr_q would have been 0);
		// - the previous slice was not conflict-free (or comb.chr_q would have been non-negative)
		for (int j = prev_slice_word_count - 1; j >= 0; --j)
		{
			bitmap_word_t AA = prev_X.AA[j] | prev_Y.AA[j] | prev_S.AA[j];
			bitmap_word_t BB = prev_X.BB[j] | prev_Y.BB[j] | prev_S.BB[j];
			bitmap_word_t conflict_mask = AA & BB;
			if (conflict_mask != 0)
			{
				// Number of consecutive no-conflict SNPs at the end of the slice
				// (see similar code above)
				prev_q = (prev_slice_word_count - 1 - j) * snp_word_bits
					+ static_cast<int>(_lzcnt_u64(conflict_mask)) - prev_slice_padding;

				break;
			}
		}
	
		// Same check yet again, but now with actual values for both prev_q and p
		if ((p < slice_snp_count || chromosome_ends) && prev_q < min_range_q[p])
			return -1;  // not calculated, starts in this slice
	}

	// Is the slice conflict-free? Increase the length of the current region,
	// which continues into the next slice (unless at the end of the chromosome)
	if (p == slice_snp_count && !chromosome_ends)
		return prev_q + p;  // region starts before this slice, and is still open

	// Not conflict-free (or at the end of a chromosome). We found the end of
	// a region that is longer than our threshold. Do something with it.
	std::vector<int> C;
	extract_final_comb(comb, C);

	std::cout << "REGION chr" << chr << ": #" << (slice_start_snp - prev_q)
		<< " to #" << (slice_start_snp + p - 1) /* inclusive */ << ", individuals ";
	const char* delim = "";
	for (int x : C)
	{
		std::cout << delim << (x + 1);
		delim = ", ";
	}
	std::cout << std::endl;

	// TODO save to file and/or store

	return -1;  // not calculated: we're either at the end of the chrosmosome so we don't care,
				// or we know a region starts in this slice, but we haven't done the work to
				// find the exact location of the last conflict in the slice
}


// Combinations
// ----------------------------------------------------------------------------

void Detector::compute_combinations()
{
	// Calculate combinations and how they are built from smaller combinations

	std::cerr << "Planning calculation of " << N << "-choose-" << K << " combinations" << std::endl;

	// Start with one-element "combinations"
	for (int i = 0; i < N; ++i)
	{
		combs.push_back({ -1, i, -1 });
		next.push_back(-1);
	}

	// Calculate combinations for all possible pairs of first, last element
	for (int x = 0; x < N - K + 1; ++x)
		for (int y = x + K - 1; y < N; ++y)
			comb(x, y, K, true);

	// Cleanup
	next.clear();
	mem.clear();

	// The above created all combinations of N choose K in final_combs,
	// and all the smaller combinations they're built from in combs.
	// The first N elements of combs are one-element "combinations", i.e.
	// the inputs themselves.

	// We allocate bitmaps only for the smaller combinations.
	bitmaps.resize(combs.size());

	std::cerr << final_combs.size() << " " << K << "-element combinations (built from "
		<< (combs.size() - N) << " smaller combinations)" << std::endl;
}


// Compute combinations of N choose K where the first element is x and the last one is y.
//
// If 'is_final' is true, store the result in 'final_combs' and assume it's
// not needed to calculate other combinations. The return value is not used.
// If 'is_final' is false, store the combinations in 'combs'. Build a
// singly-linked list of these combinations using the parallel vector 'next'
// (-1 means end of the list), and return the index of the head of this list.
// Also use 'mem' to memorize this return value for the already calculated
// (K, x, y) triples.
//
// A combination (either in 'final_combs' or in 'combs') consists of a pair
// X, Y of combinations of half the size (two indices in 'combs'). If K is
// odd, there's also a single element S (a choose-1 combination) between them
// (S is also an index in 'combs', 0 <= S < N). If K is even, S is -1.
//
int Detector::comb(int x, int y, int K, bool is_final)
{
	assert(K >= 1);
	assert(K <= y - x + 1);

	if (K == 1)
		return x;

	int tail = -1;
	int* headptr = &tail;

	if (!is_final)
	{
		// Memorize: calculate the function for the same parameters only once.
		// Insert unless the key already exists; in either case, get an iterator
		// to the new/existing item
		auto mem_result = mem.insert({ {K, x, y}, -1 });
		if (!mem_result.second)  // false means not inserted
			return mem_result.first->second;  // already calculated, we're done

		// Build a singly-linked list of combinations for this x, y and K.
		// We grow the tail rather than the head to preserve the order. It's not
		// strictly necessary, but it makes inspecting the combinations easier.
		headptr = &mem_result.first->second;
	}

	if (K == 2)
	{
		// Only one combination
		add_comb(x, -1, y, *headptr, tail, is_final);
	}
	else if (K == 3)
	{
		// First and last element are fixed; the second element can be anything between them
		for (int s = x + 1; s < y; ++s)
			add_comb(x, s, y, *headptr, tail, is_final);
	}
	else
	{
		int J = K / 2;
		if (K % 2 == 0)
		{
			// Make all combinations of two J-element halves

			// The first element of the first half must be x.
			// The last element must be such that:
			// - it has room for J elements;
			// - it leaves at least J elements for the second half.
			for (int y1 = x + J - 1; y1 < y - J + 1; ++y1)
			{
				int i0 = comb(x, y1, J, false);

				// The last element of the second half must be y.
				// The first element must be such that:
				// - it does not overlap the first half;
				// - it has room for J elements.
				for (int x1 = y1 + 1; x1 < y - J + 2; ++x1)
				{
					int j0 = comb(x1, y, J, false);

					// Each pair of first half + second half makes a K-element combination
					for (int i = i0; i >= 0; i = next[i])
						for (int j = j0; j >= 0; j = next[j])
							add_comb(i, -1, j, *headptr, tail, is_final);
				}
			}
		}
		else
		{
			// Make all combinations of two J-element halves, plus one extra between them

			// The first element of the first half must be x.
			// The last element must be such that:
			// - it has room for J elements;
			// - it leaves at least J + 1 elements for the extra and the second half.
			for (int y1 = x + J - 1; y1 < y - J; ++y1)
			{
				int i0 = comb(x, y1, J, false);

				// The last element of the second half must be y.
				// The first element must be such that:
				// - it leaves room for the extra after the first half;
				// - it has room for J elements.
				for (int x1 = y1 + 2; x1 < y - J + 2; ++x1)
				{
					int j0 = comb(x1, y, J, false);

					// Each triple of first half + extra + second half makes a K-element combination
					for (int i = i0; i >= 0; i = next[i])
						for (int s = y1 + 1; s < x1; ++s)
							for (int j = j0; j >= 0; j = next[j])
								add_comb(i, s, j, *headptr, tail, is_final);
				}
			}
		}
	}

	return *headptr;
}

void Detector::add_comb(int x, int s, int y, int& head, int& tail, bool is_final)
{
	if (is_final)
		final_combs.push_back({ x, s, y, 0 });
	else
	{
		combs.push_back({ x, s, y });
		if (tail < 0)
			tail = head = static_cast<int>(next.size());
		else
			tail = next[tail] = static_cast<int>(next.size());
		next.push_back(-1);
	}
}

void Detector::extract_comb(int comb_index, std::vector<int>& C)
{
	if (comb_index < N)
		C.push_back(comb_index);  // single element
	else
	{
		auto& comb = combs[comb_index];
		extract_comb(comb.X, C);
		if (comb.S >= 0)
			C.push_back(comb.S);
		extract_comb(comb.Y, C);
	}
}

void Detector::extract_final_comb(const final_comb_t& comb, std::vector<int>& C)
{
	extract_comb(comb.X, C);
	if (comb.S >= 0)
		C.push_back(comb.S);
	extract_comb(comb.Y, C);
}


// Command-line handling
// ----------------------------------------------------------------------------

int main(int argc, const char * argv[])
{
	std::ifstream f("t_example.tped");  // TODO command line
	Detector d;
	d.process_input(f);

	f.close();
	return 0;
}
