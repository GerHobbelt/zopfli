/*
Copyright 2011 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: lode.vandevenne@gmail.com (Lode Vandevenne)
Author: jyrki.alakuijala@gmail.com (Jyrki Alakuijala)
*/

#include "squeeze.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>

#include "blocksplitter.h"
#include "deflate.h"
#include "symbols.h"
#include "tree.h"
#include "util.h"

typedef struct SymbolStats {
  /* The literal and length symbols. */
  size_t litlens[ZOPFLI_NUM_LL];
  /* The 32 unique dist symbols, not the 32768 possible dists. */
  size_t dists[ZOPFLI_NUM_D];

  /* Length of each lit/len symbol in bits. */
  double ll_symbols[ZOPFLI_NUM_LL];
  /* Length of each dist symbol in bits. */
  double d_symbols[ZOPFLI_NUM_D];
} SymbolStats;

/* Sets everything to 0. */
static void InitStats(SymbolStats* stats) {
  memset(stats->litlens, 0, ZOPFLI_NUM_LL * sizeof(stats->litlens[0]));
  memset(stats->dists, 0, ZOPFLI_NUM_D * sizeof(stats->dists[0]));

  memset(stats->ll_symbols, 0, ZOPFLI_NUM_LL * sizeof(stats->ll_symbols[0]));
  memset(stats->d_symbols, 0, ZOPFLI_NUM_D * sizeof(stats->d_symbols[0]));
}

static void CopyStats(SymbolStats* source, SymbolStats* dest) {
  memcpy(dest->litlens, source->litlens,
         ZOPFLI_NUM_LL * sizeof(dest->litlens[0]));
  memcpy(dest->dists, source->dists, ZOPFLI_NUM_D * sizeof(dest->dists[0]));

  memcpy(dest->ll_symbols, source->ll_symbols,
         ZOPFLI_NUM_LL * sizeof(dest->ll_symbols[0]));
  memcpy(dest->d_symbols, source->d_symbols,
         ZOPFLI_NUM_D * sizeof(dest->d_symbols[0]));
}

/* Adds the bit lengths. */
static void AddWeighedStatFreqs(const SymbolStats* stats1, double w1,
                                const SymbolStats* stats2, double w2,
                                SymbolStats* result) {
  size_t i;
  for (i = 0; i < ZOPFLI_NUM_LL; i++) {
    result->litlens[i] =
        (size_t) (stats1->litlens[i] * w1 + stats2->litlens[i] * w2);
  }
  for (i = 0; i < ZOPFLI_NUM_D; i++) {
    result->dists[i] =
        (size_t) (stats1->dists[i] * w1 + stats2->dists[i] * w2);
  }
  result->litlens[256] = 1;  /* End symbol. */
}

typedef struct RanState {
  unsigned int m_w, m_z;
} RanState;

static void InitRanState(RanState* state) {
  state->m_w = 1;
  state->m_z = 2;
}

/* Get random number: "Multiply-With-Carry" generator of G. Marsaglia */
static unsigned int Ran(RanState* state) {
  state->m_z = 36969 * (state->m_z & 65535) + (state->m_z >> 16);
  state->m_w = 18000 * (state->m_w & 65535) + (state->m_w >> 16);
  return (state->m_z << 16) + state->m_w;  /* 32-bit result. */
}

static void RandomizeFreqs(RanState* state, size_t* freqs, int n, int every, int out_of) {
  int i;
  for (i = 0; i < n; i++) {
    if ((Ran(state) >> 4) % out_of < every) freqs[i] = freqs[Ran(state) % n];
  }
}

/* Switch *every* N out of *out_of* elements*/
static void RandomizeStatFreqs(RanState* state, SymbolStats* stats, int every, int out_of) {
  RandomizeFreqs(state, stats->litlens, ZOPFLI_NUM_LL, every, out_of);
  RandomizeFreqs(state, stats->dists, ZOPFLI_NUM_D, every, out_of);
  stats->litlens[256] = 1;  /* End symbol. */
}

static void ClearStatFreqs(SymbolStats* stats) {
  size_t i;
  for (i = 0; i < ZOPFLI_NUM_LL; i++) stats->litlens[i] = 0;
  for (i = 0; i < ZOPFLI_NUM_D; i++) stats->dists[i] = 0;
}

/*
Function that calculates a cost based on a model for the given LZ77 symbol.
litlen: means literal symbol if dist is 0, length otherwise.
*/
typedef double CostModelFun(unsigned litlen, unsigned dist, void* context);

/*
Cost model which should exactly match fixed tree.
type: CostModelFun
*/
static double GetCostFixed(unsigned litlen, unsigned dist, void* unused) {
  (void)unused;
  if (dist == 0) {
    if (litlen <= 143) return 8;
    else return 9;
  } else {
    int dbits = ZopfliGetDistExtraBits(dist);
    int lbits = ZopfliGetLengthExtraBits(litlen);
    int lsym = ZopfliGetLengthSymbol(litlen);
    int cost = 0;
    if (lsym <= 279) cost += 7;
    else cost += 8;
    cost += 5;  /* Every dist symbol has length 5. */
    return cost + dbits + lbits;
  }
}

/*
Cost model based on symbol statistics.
type: CostModelFun
*/
static double GetCostStat(unsigned litlen, unsigned dist, void* context) {
  SymbolStats* stats = (SymbolStats*)context;
  if (dist == 0) {
    return stats->ll_symbols[litlen];
  } else {
    int lsym = ZopfliGetLengthSymbol(litlen);
    int lbits = ZopfliGetLengthExtraBits(litlen);
    int dsym = ZopfliGetDistSymbol(dist);
    int dbits = ZopfliGetDistExtraBits(dist);
    return lbits + dbits + stats->ll_symbols[lsym] + stats->d_symbols[dsym];
  }
}

/*
Finds the minimum possible cost this cost model can return for valid length and
distance symbols.
*/
static double GetCostModelMinCost(CostModelFun* costmodel, void* costcontext) {
  double mincost;
  int bestlength = 0; /* length that has lowest cost in the cost model */
  int bestdist = 0; /* distance that has lowest cost in the cost model */
  int i;
  /*
  Table of distances that have a different distance symbol in the deflate
  specification. Each value is the first distance that has a new symbol. Only
  different symbols affect the cost model so only these need to be checked.
  See RFC 1951 section 3.2.5. Compressed blocks (length and distance codes).
  */
  static const int dsymbols[30] = {
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513,
    769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
  };

  mincost = ZOPFLI_LARGE_FLOAT;
  for (i = 3; i < 259; i++) {
    double c = costmodel(i, 1, costcontext);
    if (c < mincost) {
      bestlength = i;
      mincost = c;
    }
  }

  mincost = ZOPFLI_LARGE_FLOAT;
  for (i = 0; i < 30; i++) {
    double c = costmodel(3, dsymbols[i], costcontext);
    if (c < mincost) {
      bestdist = dsymbols[i];
      mincost = c;
    }
  }

  return costmodel(bestlength, bestdist, costcontext);
}

static size_t zopfli_min(size_t a, size_t b) {
  return a < b ? a : b;
}

static double zopfli_dmin(double a, double b) {
  return a < b ? a : b;
}

/*
Performs the forward pass for "squeeze". Gets the most optimal length to reach
every byte from a previous byte, using cost calculations.
s: the ZopfliBlockState
in: the input data array
instart: where to start
inend: where to stop (not inclusive)
costmodel: function to calculate the cost of some lit/len/dist pair.
costcontext: abstract context for the costmodel function
length_array: output array of size (inend - instart) which will receive the best
    length to reach this byte from a previous byte.
returns the cost that was, according to the costmodel, needed to get to the end.
*/
static double GetBestLengths(ZopfliBlockState *s,
                             const unsigned char* in,
                             size_t instart, size_t inend,
                             CostModelFun* costmodel, void* costcontext,
                             unsigned short* length_array,
                             ZopfliHash* h, float* costs) {
  /* Best cost to get here so far. */
  size_t blocksize = inend - instart;
  size_t i = 0, k, kend;
  unsigned short leng;
  unsigned short dist;
  unsigned short sublen[259];
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;
  double result;
  double mincost = GetCostModelMinCost(costmodel, costcontext);
  double mincostaddcostj;

  if (instart == inend) return 0;

  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) {
    ZopfliUpdateHash(in, i, inend, h);
  }

  for (i = 1; i < blocksize + 1; i++) costs[i] = ZOPFLI_LARGE_FLOAT;
  costs[0] = 0;  /* Because it's the start. */
  length_array[0] = 0;

  for (i = instart; i < inend; i++) {
    size_t j = i - instart;  /* Index in the costs array and length_array. */
    ZopfliUpdateHash(in, i, inend, h);

#ifdef ZOPFLI_SHORTCUT_LONG_REPETITIONS
    /* If we're in a long repetition of the same character and have more than
    ZOPFLI_MAX_MATCH characters before and after our position. */
    if (h->same[i & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH * 2
        && i > instart + ZOPFLI_MAX_MATCH + 1
        && i + ZOPFLI_MAX_MATCH * 2 + 1 < inend
        && h->same[(i - ZOPFLI_MAX_MATCH) & ZOPFLI_WINDOW_MASK]
            > ZOPFLI_MAX_MATCH) {
      double symbolcost = costmodel(ZOPFLI_MAX_MATCH, 1, costcontext);
      /* Set the length to reach each one to ZOPFLI_MAX_MATCH, and the cost to
      the cost corresponding to that length. Doing this, we skip
      ZOPFLI_MAX_MATCH values to avoid calling ZopfliFindLongestMatch. */
      for (k = 0; k < ZOPFLI_MAX_MATCH; k++) {
        costs[j + ZOPFLI_MAX_MATCH] = costs[j] + symbolcost;
        length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH;
        i++;
        j++;
        ZopfliUpdateHash(in, i, inend, h);
      }
    }
#endif

    ZopfliFindLongestMatch(s, h, in, i, inend, ZOPFLI_MAX_MATCH, sublen,
                           &dist, &leng);

    /* Literal. */
    if (i + 1 <= inend) {
      double newCost = costmodel(in[i], 0, costcontext) + costs[j];
      assert(newCost >= 0);
      if (newCost < costs[j + 1]) {
        costs[j + 1] = newCost;
        length_array[j + 1] = 1;
      }
    }
    /* Lengths. */
    kend = zopfli_min(leng, inend-i);
    mincostaddcostj = mincost + costs[j];
    for (k = 3; k <= kend; k++) {
      double newCost;

      /* Calling the cost model is expensive, avoid this if we are already at
      the minimum possible cost that it can return. */
     if (costs[j + k] <= mincostaddcostj) continue;

      newCost = costmodel(k, sublen[k], costcontext) + costs[j];
      assert(newCost >= 0);
      if (newCost < costs[j + k]) {
        assert(k <= ZOPFLI_MAX_MATCH);
        costs[j + k] = newCost;
        length_array[j + k] = k;
      }
    }
  }

  assert(costs[blocksize] >= 0);
  result = costs[blocksize];

  return result;
}

/*
Calculates the optimal path of lz77 lengths to use, from the calculated
length_array. The length_array must contain the optimal length to reach that
byte. The path will be filled with the lengths to use, so its data size will be
the amount of lz77 symbols.
*/
static void TraceBackwards(size_t size, const unsigned short* length_array,
                           unsigned short** path, size_t* pathsize) {
  size_t index = size;
  if (size == 0) return;
  for (;;) {
    ZOPFLI_APPEND_DATA(length_array[index], path, pathsize, ZOPFLI_DYN_ALLOC);
    assert(length_array[index] <= index);
    assert(length_array[index] <= ZOPFLI_MAX_MATCH);
    assert(length_array[index] != 0);
    index -= length_array[index];
    if (index == 0) break;
  }

  /* Mirror result. */
  for (index = 0; index < *pathsize / 2; index++) {
    unsigned short temp = (*path)[index];
    (*path)[index] = (*path)[*pathsize - index - 1];
    (*path)[*pathsize - index - 1] = temp;
  }
}

static void FollowPath(ZopfliBlockState* s,
                       const unsigned char* in, size_t instart, size_t inend,
                       unsigned short* path, size_t pathsize,
                       ZopfliLZ77Store* store, ZopfliHash *h) {
  size_t i, j, pos = 0;
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;

  size_t total_length_test = 0;

  if (instart == inend) return;

  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) {
    ZopfliUpdateHash(in, i, inend, h);
  }

  pos = instart;
  for (i = 0; i < pathsize; i++) {
    unsigned short length = path[i];
    unsigned short dummy_length;
    unsigned short dist;
    assert(pos < inend);

    ZopfliUpdateHash(in, pos, inend, h);

    /* Add to output. */
    if (length >= ZOPFLI_MIN_MATCH) {
      /* Get the distance by recalculating longest match. The found length
      should match the length from the path. */
      ZopfliFindLongestMatch(s, h, in, pos, inend, length, 0,
                             &dist, &dummy_length);
      assert(!(dummy_length != length && length > 2 && dummy_length > 2));
      ZopfliVerifyLenDist(in, inend, pos, dist, length);
      ZopfliStoreLitLenDist(length, dist, pos, store);
      total_length_test += length;
    } else {
      length = 1;
      ZopfliStoreLitLenDist(in[pos], 0, pos, store);
      total_length_test++;
    }


    assert(pos + length <= inend);
    for (j = 1; j < length; j++) {
      ZopfliUpdateHash(in, pos + j, inend, h);
    }

    pos += length;
  }
}

/* Calculates the entropy of the statistics */
static void CalculateStatistics(SymbolStats* stats) {
  ZopfliCalculateEntropy(stats->litlens, ZOPFLI_NUM_LL, stats->ll_symbols);
  ZopfliCalculateEntropy(stats->dists, ZOPFLI_NUM_D, stats->d_symbols);
}

/* Appends the symbol statistics from the store. */
static void GetStatistics(const ZopfliLZ77Store* store, SymbolStats* stats) {
  size_t i;
  for (i = 0; i < store->size; i++) {
    if (store->dists[i] == 0) {
      stats->litlens[store->litlens[i]]++;
    } else {
      stats->litlens[ZopfliGetLengthSymbol(store->litlens[i])]++;
      stats->dists[ZopfliGetDistSymbol(store->dists[i])]++;
    }
  }
  stats->litlens[256] = 1;  /* End symbol. */

  CalculateStatistics(stats);
}

/*
Does a single run for ZopfliLZ77Optimal. For good compression, repeated runs
with updated statistics should be performed.
s: the block state
in: the input data array
instart: where to start
inend: where to stop (not inclusive)
path: pointer to dynamically allocated memory to store the path
pathsize: pointer to the size of the dynamic path array
length_array: array of size (inend - instart) used to store lengths
costmodel: function to use as the cost model for this squeeze run
costcontext: abstract context for the costmodel function
store: place to output the LZ77 data
returns the cost that was, according to the costmodel, needed to get to the end.
    This is not the actual cost.
*/
static double LZ77OptimalRun(ZopfliBlockState* s,
    const unsigned char* in, size_t instart, size_t inend,
    unsigned short** path, size_t* pathsize,
    unsigned short* length_array, CostModelFun* costmodel,
    void* costcontext, ZopfliLZ77Store* store,
    ZopfliHash* h, float* costs) {

  double cost = GetBestLengths(s, in, instart, inend, costmodel,
                costcontext, length_array, h, costs);
  free(*path);
  *path = 0;
  *pathsize = 0;
  TraceBackwards(inend - instart, length_array, path, pathsize);
  FollowPath(s, in, instart, inend, *path, *pathsize, store, h);
  assert(cost < ZOPFLI_LARGE_FLOAT);
  return cost;
}

/*
 * OMP Parallel parameters
 */
#ifndef AM_OMP_THREAD_NUM
  #define AM_OMP_THREAD_NUM 4
#endif

#ifndef AM_OMP_T_RAND_NUM
  #define AM_OMP_T_RAND_NUM 1 
#endif

#ifndef AM_OMP_T_RAND_DENOM
  #define AM_OMP_T_RAND_DENOM 20
#endif

void ZopfliLZ77Optimal(ZopfliBlockState *s,
                       const unsigned char* in, size_t instart, size_t inend,
                       int numiterations,
                       ZopfliLZ77Store* store) {
  /* Dist to get to here with smallest cost. */
  const size_t blocksize = inend - instart;
  
  SymbolStats beststats;
  int i=0;
  double bestcost = ZOPFLI_LARGE_FLOAT;
  double best_iter_cost= ZOPFLI_LARGE_FLOAT; /* best cost of the iteration */
  SymbolStats * broadcast_stats;
  /*
   * thread private declarations initialized up here for C90 reasons
   */
  int t_i;
  ZopfliLZ77Store t_currentstore;
  ZopfliHash t_hash;
  SymbolStats t_stats,  t_laststats;
  unsigned short* t_path;
  size_t t_pathsize;
  double t_cost;
  double t_lastcost;
  int t_my_tid;
  unsigned short* t_length_array;
  RanState t_ran_state;

  int t_lastrandomstep;
  float* t_costs;
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  assert(0); /* not allowed with OMP */
#endif 

  /* TODO s is the block state -- should this be parallel or copied? */

  /* Allocate thread private variables*/
  #pragma omp parallel num_threads(AM_OMP_THREAD_NUM) \
    default(none) \
    private(t_length_array, t_path, t_pathsize, t_currentstore, t_hash, \
      t_stats, t_laststats, t_cost, t_lastcost, t_ran_state, t_lastrandomstep, \
      t_costs, t_my_tid, t_i) \
    shared(beststats, i, bestcost, broadcast_stats, best_iter_cost, in, instart, \
            inend, numiterations, s, stderr, store)
  {
    t_length_array = (unsigned short*)malloc(sizeof(unsigned short) * (blocksize + 1));
    t_path = 0;
    t_pathsize = 0;
    t_cost = ZOPFLI_LARGE_FLOAT;
    t_lastcost = ZOPFLI_LARGE_FLOAT;
    t_my_tid = omp_get_thread_num();

    /* Try randomizing the costs a bit off the bat to take advantage of the parallelism */
    t_lastrandomstep = -1;
    t_costs = (float*)malloc(sizeof(float) * (blocksize + 1));

    if (!t_costs) exit(-1); /* Allocation failed. */
    if (!t_length_array) exit(-1); /* Allocation failed. */

    InitRanState(&t_ran_state);
    /* TODO Add TID to each ran_state variable for divergence*/
    InitStats(&t_stats);
    ZopfliInitLZ77Store(in, &t_currentstore);
    ZopfliAllocHash(ZOPFLI_WINDOW_SIZE, &t_hash);

    /* Do regular deflate on thread 0, 
    then loop multiple batches of shortest path runs with randomized stats,
    each time using the best statistics of the previous batch. */

    /* Initial run. */
    if (t_my_tid == 0) {
      ZopfliLZ77Greedy(s, in, instart, inend, &t_currentstore, &t_hash);
      GetStatistics(&t_currentstore, &t_stats);
      memcpy(&beststats, &t_stats, sizeof(t_stats));
      /* broadcast stats*/
      broadcast_stats = &t_stats;
    }
    #pragma omp barrier
      /* receive stats -- is there a better way to do this with a copy clause for a struct? */
    #pragma omp parallel for
    for (i = 0; i < AM_OMP_THREAD_NUM; i++) {
      if (&t_stats != broadcast_stats) {
        memcpy(&t_stats, broadcast_stats, sizeof(t_stats));
      }
    }

    #pragma omp barrier

    /* Repeat statistics with each time the cost model from the previous stat run. 
      Each thread runs t_i times TODO - AM drop the number of iterations*/
    for (t_i = 0; t_i < numiterations; t_i++) {

      /* Randomize the statistics slightly for each thread */
      #pragma omp parallel for
      for(i = 0; i < AM_OMP_THREAD_NUM; i++) {
        CopyStats(&beststats, &t_stats);
        /* Randomization is dependent on the thread - try 1 and 20*/
        RandomizeStatFreqs(&t_ran_state, &t_stats, t_my_tid * AM_OMP_T_RAND_NUM, AM_OMP_T_RAND_DENOM);
        CalculateStatistics(&t_stats);
        t_lastrandomstep = i;
      }

      /* Run the Optimal run with slightly different stats on each thread. */
      #pragma omp parallel for
      for(i = 0; i < AM_OMP_THREAD_NUM; i++) {
        DEBUG_PRINT(("starting optimal course thread,ti,i:%d,%d,%d \n", t_my_tid, t_i, i));
        ZopfliCleanLZ77Store(&t_currentstore);
        ZopfliInitLZ77Store(in, &t_currentstore);
        LZ77OptimalRun(s, in, instart, inend, &t_path, &t_pathsize,
                    t_length_array, GetCostStat, (void*)&t_stats,
                    &t_currentstore, &t_hash, t_costs);
        t_cost = ZopfliCalculateBlockSize(&t_currentstore, 0, t_currentstore.size, 2);
        DEBUG_PRINT(("Calculated Cost: threads,ti,i,cost: %d,%d,%d,%f\n", t_my_tid, t_i, i, t_cost));
        if (s->options->verbose_more || (s->options->verbose && t_cost < bestcost)) {
            (void)0;
            DEBUG_PRINT(("Found new best cost: threads,ti,i,cost: %d,%d,%d,%f\n", t_my_tid, t_i, i, t_cost));
            /*
            fprintf(stderr, "Calculated Cost: threads,ti,i,cost: %d,%d %d %f\n", t_my_tid, t_i, i, t_cost);
            fprintf(stderr, "threads, %d, tpath %p, tpathsize %lu, tl\n", t_my_tid, t_path, t_pathsize);
            */
        }
      } /* OMP barrier*/

       best_iter_cost= ZOPFLI_LARGE_FLOAT; /* best cost of the iteration */
      #pragma omp barrier
      /* OMP Min Reduce*/
      DEBUG_PRINT(("About to reduce: thread,ti,cost: %d,%d,%f\n", t_my_tid, t_i, t_cost));
      #pragma omp parallel for reduction(min:best_iter_cost)
      for(i = 0; i < AM_OMP_THREAD_NUM; i++) {
          /* TODO this is not reducing (or if it is, its reducing too well) */
        best_iter_cost = zopfli_dmin(best_iter_cost, t_cost);
      }

      /* Format spacing */
      DEBUG_PRINT(("\n"));

      /* OMP Min Reduce*/
      #pragma omp parallel
      {
      DEBUG_PRINT(("evaluating cost to bestitercost(%f): thread,Iteration,cost:%d,%d,%f\n",best_iter_cost, t_my_tid, t_i, t_cost));
        if (t_cost == best_iter_cost && t_cost < bestcost) {
            DEBUG_PRINT(("Entering CriticalZone: thread,Iteration,cost:%d,%d,%f\n", t_my_tid, t_i, t_cost));
            #pragma omp critical (copy_to_best)
            {
                /* TODO Can this be merged with the last loop?*/
                /* Copy to the output store. */
                fprintf(stderr, "threads,Iteration: %d,ti %d copying new best cost %f\n", t_my_tid, t_i, t_cost);
                ZopfliCopyLZ77Store(&t_currentstore, store);
                /** This! don't do this */
                CopyStats(&t_stats, &beststats);
                bestcost = t_cost;
            }
            DEBUG_PRINT(("Leaving CriticalZone: thread,Iteration,cost:%d,%d,%f\n", t_my_tid, t_i, t_cost));
        }
      }

      CopyStats(&t_stats, &t_laststats);
      ClearStatFreqs(&t_stats);
      GetStatistics(&t_currentstore, &t_stats);

      /* Format spacing */
      DEBUG_PRINT(("\n"));

      /* TODO experiment with randomizing and doing the last weighted stats variable.*/
      if (t_lastrandomstep != -1) {
        /* This makes it converge slower but better. Do it only once the
        randomness kicks in so that if the user does few iterations, it gives a
        better result sooner. */
        AddWeighedStatFreqs(&t_stats, 1.0, &t_laststats, 0.5, &t_stats);
        CalculateStatistics(&t_stats);
      }

      /* Help figure out the convergance effectiveness.*/
      DEBUG_PRINT(("thread,Iteration:%d,%d, cost: %f, bestCost: %f\n", t_my_tid, i, t_cost, bestcost));
      t_lastcost = t_cost;
    }

  free(t_length_array);
  free(t_path);
  free(t_costs);
  ZopfliCleanLZ77Store(&t_currentstore);
  ZopfliCleanHash(&t_hash);
  } /* omp parallel*/

  /* TODO print stats to get a frequency distribution ish for montecarlo? */
#ifdef DEBUG
  DEBUG_PRINT((" litlens:\n"));
  for(i = 0; i < ZOPFLI_NUM_LL; i++) {
      DEBUG_PRINT(("%zu ", beststats.litlens[i]));
  }
  DEBUG_PRINT(("\ndists:\n"));
  for(i = 0; i < ZOPFLI_NUM_D; i++) {
      DEBUG_PRINT(("%zu ", beststats.dists[i]));
  }
  DEBUG_PRINT(("\n\n"));
#endif
}

void ZopfliLZ77OptimalFixed(ZopfliBlockState *s,
                            const unsigned char* in,
                            size_t instart, size_t inend,
                            ZopfliLZ77Store* store)
{
  /* Dist to get to here with smallest cost. */
  size_t blocksize = inend - instart;
  unsigned short* length_array =
      (unsigned short*)malloc(sizeof(unsigned short) * (blocksize + 1));
  unsigned short* path = 0;
  size_t pathsize = 0;
  ZopfliHash hash;
  ZopfliHash* h = &hash;
  float* costs = (float*)malloc(sizeof(float) * (blocksize + 1));

  if (!costs) exit(-1); /* Allocation failed. */
  if (!length_array) exit(-1); /* Allocation failed. */

  ZopfliAllocHash(ZOPFLI_WINDOW_SIZE, h);

  s->blockstart = instart;
  s->blockend = inend;

  /* Shortest path for fixed tree This one should give the shortest possible
  result for fixed tree, no repeated runs are needed since the tree is known. */
  LZ77OptimalRun(s, in, instart, inend, &path, &pathsize,
                 length_array, GetCostFixed, 0, store, h, costs);

  free(length_array);
  free(path);
  free(costs);
  ZopfliCleanHash(h);

}
