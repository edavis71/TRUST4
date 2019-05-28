// The data structure holds the set of sequences (can be "assembled" from several reads)
#ifndef _MOURISL_SEQSET_HEADER
#define _MOURISL_SEQSET_HEADER

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <queue>

#include "SimpleVector.hpp"
#include "KmerIndex.hpp"
#include "ReadFiles.hpp"
#include "AlignAlgo.hpp"


struct _seqWrapper
{
	char *name ;
	char *consensus ; // This should be handled by malloc/free.
	int consensusLen ;
	SimpleVector<struct _posWeight> posWeight ;
	bool isRef ; // this is from reference.

	int minLeftExtAnchor, minRightExtAnchor ; // only overlap with size larger than this can be counted as valid extension.

	struct _pair info[3] ; // For storing extra information. for ref, info[0,1] contains the coordinate for CDR1
} ;

struct _hit
{
	struct _indexInfo indexHit ;
	int readOffset ;
	int strand ; // -1: different strand, 1: same strand. When strand==-1, the readOffset is the offset in the rcSeq.
	
	bool operator<( const struct _hit &b ) const
	{
		if ( strand != b.strand )
			return strand < b.strand ;
		else if ( indexHit.idx != b.indexHit.idx )
			return indexHit.idx < b.indexHit.idx ;
		else if ( readOffset != b.readOffset )
			return  readOffset < b.readOffset ;
		else if ( indexHit.offset != b.indexHit.offset )
			return indexHit.offset < b.indexHit.offset ;
		
		return false ;
	}
} ;

struct _overlap
{
	int seqIdx ;
	int readStart, readEnd ; // When strand ==-1, the start,end is the offset in the rcSeq.
	int seqStart, seqEnd ;
	int strand ;
	
	int matchCnt ; // The number of matched bases, count TWICE.
	double similarity ;

	SimpleVector<struct _pair> *hitCoords ;
	SimpleVector<int> *info ; // store extra informations 

	bool operator<( const struct _overlap &b ) const
	{
		// The overlap with more matched bases should come first.
		if ( matchCnt > b.matchCnt + 2 || matchCnt < b.matchCnt - 2 )
			return matchCnt > b.matchCnt ;
		else if ( similarity != b.similarity )
			return similarity > b.similarity ; 
		else if ( readEnd - readStart != b.readEnd - b.readStart )
			return readEnd - readStart > b.readEnd - b.readStart ;
		else if ( seqIdx != b.seqIdx )
			return seqIdx < b.seqIdx ;
		else if ( strand !=  b.strand )
			return strand < b.strand ;
		else if ( readStart != b.readStart )
			return readStart < b.readStart ;
		else if ( readEnd != b.readEnd )
			return readEnd < b.readEnd ;
		else if ( seqStart != b.seqStart )
			return seqStart < b.seqStart ;
		else 
			return seqEnd < b.seqEnd ; 

		return false ;
	}
} ;

struct _assignRead
{
	char *id ;
	char *read ;

	struct _overlap overlap ;
} ;


class SeqSet
{
private:
	std::vector<struct _seqWrapper> seqs ;
	KmerIndex seqIndex ;
	int kmerLength ;
	int radius ;
	int hitLenRequired ;

	// Some threshold
	double novelSeqSimilarity ;
	double refSeqSimilarity ;
	double repeatSimilarity ; // e.g., the repeat when building the branch graph.

	struct _overlap prevAddInfo ; 

	static bool CompSortPairBInc( const struct _pair &p1, const struct _pair &p2 )
	{
		if ( p1.b != p2.b )
			return p1.b < p2.b ;
		else
			return p1.a < p2.a ;
	}
	
	static bool CompSortPairAInc( const struct _pair &p1, const struct _pair &p2 )
	{
		return p1.a < p2.a ;
	}

	static bool CompSortOverlapsOnReadCoord( const struct _overlap &a, const struct _overlap &b )
	{
		return a.readStart < b.readStart ; 
	}

	static bool CompSortAssignedReadById( const struct _assignRead &a, const struct _assignRead &b )
	{
		return strcmp( a.id, b.id ) < 0 ;
	}
		
	
	bool IsReverseComplement( char *a, char *b )
	{
		int i, j ;
		int len = strlen( a ) ;
		if ( len != strlen( b) ) 
			return false ;
		for ( i = 0, j = len - 1 ; i < len ; ++i, --j )
			if ( a[i] == 'N' && b[j] == 'N' )
				continue ;
			else if ( a[i] != 'N' && b[j] != 'N' )
			{
				if ( 3 - nucToNum[ a[i] - 'A' ] != nucToNum[ b[j] - 'A' ] )
					return false ;
			}
			else
				return false ;
		return true ;
	}

	void Reverse( char *r, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
			r[i] = seq[len - 1 - i] ;  
		r[i] = '\0' ;
	}
	
	bool IsPosWeightCompatible( const struct _posWeight &a, const struct _posWeight &b )
	{
		int sumA = a.Sum() ;
		int sumB = b.Sum() ;
	
		if ( sumA == 0 || sumB == 0 
			|| ( sumA < 3 * a.count[0] && sumB < 3 * b.count[0] ) 
			|| ( sumA < 3 * a.count[1] && sumB < 3 * b.count[1] )
			|| ( sumA < 3 * a.count[2] && sumB < 3 * b.count[2] )
			|| ( sumA < 3 * a.count[3] && sumB < 3 * b.count[3] ) )
			return true ;
		return false ;	
	}

	bool IsOverlapIntersect( const struct _overlap &a, const struct _overlap &b )
	{
		if ( a.seqIdx == b.seqIdx && 
			( ( a.seqStart <= b.seqStart && a.seqEnd >= b.seqStart ) 
			|| ( b.seqStart <= a.seqStart && b.seqEnd >= a.seqStart ) ) )
			return true ;
		return false ;
	}

	// Return the first index whose hits.a is smaller or equal to valA
	int BinarySearch_LIS( int top[], int size, int valA, SimpleVector<struct _pair> &hits )
	{
		int l = 0, r = size - 1 ;
		int m ;
		while ( l <= r )
		{
			m = ( l + r ) / 2 ;
			if ( valA == hits[ top[m] ].a )
			{
				return m ;
			}
			else if ( valA < hits[ top[m] ].a )
			{
				r = m - 1 ;
			}
			else
			{
				l = m + 1 ;
			}
		}
		return l - 1 ;
	}

	// The O(nlogn) method for solving LIS problem, suppose there are n hits.
	// Return the LIS, the LIS's length is returned by the function
	int LongestIncreasingSubsequence( SimpleVector<struct _pair> &hits, SimpleVector<struct _pair> &LIS ) 
	{
		// Only use the first hit of each qhit
		// Bias towards left

		int i, j, k ;
		int ret = 0 ;
		int size = hits.Size() ;

		int *record = new int[size] ; // The index of the selected hits
		int *top = new int[size] ; // record the index of the hits with smallest valB of the corresponding LIS length. kind of the top element.
		int *link = new int[size] ; // used to retrieve the LIS

		int rcnt = 1 ;
		record[0] = 0 ;
		for ( i = 1 ; i < size ; ++i )
		{
			//if ( hits[i].b == hits[i - 1].b )
			//	continue ;
			record[rcnt] = i ;
			++rcnt ;
		}
		top[0] = 0 ;
		link[0] = -1 ;
		ret = 1 ;
		for ( i = 1 ; i < rcnt ; ++i )
		{
			int tag = 0 ;
			if ( hits[ top[ ret - 1 ] ].a <= hits[ record[i] ].a )
				tag = ret - 1 ;
			else
				tag = BinarySearch_LIS( top, ret, hits[ record[i] ].a, hits ) ;			
			
			if ( tag == -1 )
			{
				top[0] = record[i] ;
				link[ record[i] ] = -1 ;
			}
			else if ( hits[ record[i] ].a > hits[ top[tag] ].a )
			{
				if ( tag == ret - 1 )
				{
					top[ret] = record[i] ;
					++ret ;
					link[ record[i] ] = top[tag] ;
				}
				else if ( hits[ record[i] ].a < hits[ top[tag + 1] ].a )
				{
					top[ tag + 1 ] = record[i] ;
					link[ record[i] ] = top[tag] ;
				}
			}
		}


		k = top[ret - 1] ;
		for ( i = ret - 1 ; i >= 0 ; --i )
		{
			LIS.PushBack( hits[k] ) ;
			k = link[k] ;	
		}
		LIS.Reverse() ;
		//for ( i = 0 ; i < ret ; ++i )
		//	LIS.PushBack( hits[ top[i] ] ) ;
		
		// Remove elements with same b.
		if ( ret > 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < ret ; ++i )
			{
				if ( LIS[i].b == LIS[k - 1].b )
					continue ;
				LIS[k] = LIS[i] ;
				++k ;
			}
			ret = k ;
		}

		delete []top ;
		delete []record ;
		delete []link ;

		return ret ;
	}

	void GetAlignStats( char *align, bool update, int &matchCnt, int &mismatchCnt, int &indelCnt)
	{
		int k ;
		if ( !update )
		{
			matchCnt = mismatchCnt = indelCnt = 0 ;
		}

		for ( k = 0 ; align[k] != -1 ; ++k )
		{
			if ( align[k] == EDIT_MATCH )
				++matchCnt ;
			else if ( align[k] == EDIT_MISMATCH )
				++mismatchCnt ;
			else 
				++indelCnt ;
		}
	}


	bool IsOverlapLowComplex( char *r, struct _overlap &o )
	{
		int cnt[4] = {0, 0, 0, 0} ;
		int i ;
		for ( i = o.readStart ; i <= o.readEnd ; ++i )
		{
			if ( r[i] == 'N' )
				continue ;
			++cnt[ nucToNum[ r[i] - 'A' ] ] ;
		}
		int len = o.readEnd - o.readStart + 1 ;
		int lowCnt = 0 ; 
		for ( i = 0 ; i < 4 ; ++i )
			if ( cnt[i] <= 2 )
				++lowCnt ;
		if ( lowCnt >= 2 )
			return true ;
		return false ;
	}

	void SetPrevAddInfo( int seqIdx, int readStart, int readEnd, int seqStart, int seqEnd, int strand )
	{
		prevAddInfo.seqIdx = seqIdx ;
		prevAddInfo.readStart = readStart ;
		prevAddInfo.readEnd = readEnd ;
		prevAddInfo.seqStart = seqStart ;
		prevAddInfo.seqEnd = seqEnd ;
		prevAddInfo.strand = strand ;
	}
	
	char DnaToAa( char a, char b, char c )
	{
		if ( a == 'N' || b == 'N' || c == 'N' )
			return '-' ;

		if ( a == 'A' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'K' ;
				else
					return 'N' ;
			}
			else if ( b == 'C' )
			{
				return 'T' ;
			}
			else if ( b == 'G' )
			{
				if ( c == 'A' || c == 'G' )
					return 'R' ;
				else
					return 'S' ;
			}
			else
			{
				if ( c == 'G' )
					return 'M' ;
				else 
					return 'I' ;
			}
		}
		else if ( a == 'C' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'Q' ;
				else
					return 'H' ;
				
			}
			else if ( b == 'C' )
			{
				return 'P' ;
			}
			else if ( b == 'G' )
			{
				return 'R' ;
			}
			else
			{
				return 'L' ;
			}
		}
		else if ( a == 'G' )
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return 'E' ;
				else
					return 'D' ;
			}
			else if ( b == 'C' )
			{
				return 'A' ;
			}
			else if ( b == 'G' )
			{
				return 'G' ;
			}
			else
			{
				return 'V' ;
			}
		}
		else
		{
			if ( b == 'A' )
			{
				if ( c == 'A' || c == 'G' )
					return '*' ;
				else
					return 'Y' ;
			}
			else if ( b == 'C' )
			{
				return 'S' ;
			}
			else if ( b == 'G' )
			{
				if ( c == 'A' )
					return '*' ;
				else if ( c == 'G' )
					return 'W' ;
				else
					return 'C' ;
				
			}
			else
			{
				if ( c == 'A' || c == 'G' )
					return 'L' ;
				else
					return 'F' ;
			}
		}
	}

	void ReleaseSeq( int idx )
	{
		if ( seqs[idx].consensus == NULL )
			return ;
		free( seqs[idx].name ) ;
		free( seqs[idx].consensus ) ;
		seqs[idx].posWeight.Release() ;

		seqs[idx].name = seqs[idx].consensus = NULL ;
	}
	
	// Use the hits to extract overlaps from SeqSet
	int GetOverlapsFromHits( SimpleVector<struct _hit> &hits, int hitLenRequired, int filter, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		int hitSize = hits.Size() ;
		
		SimpleVector<struct _pair> hitCoordDiff ;
		hitCoordDiff.Reserve( hitSize ) ;
		SimpleVector<struct _pair> concordantHitCoord ;
		SimpleVector<struct _pair> hitCoordLIS ;
		SimpleVector<struct _hit> finalHits ;
		
		// Compute the minHitRequired.
		int novelMinHitRequired = 3 ;
		int refMinHitRequired = 3 ;
		if ( filter == 1 )
		{
			int possibleOverlapCnt = 0 ;
			int longestHits = 0 ;
			for ( i = 0 ; i < hitSize ; ++i )
			{
				for ( j = i + 1 ; j < hitSize ; ++j )
					if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
						break ;
				if ( !seqs[ hits[i].indexHit.idx].isRef )
				{
					if ( j - i > novelMinHitRequired )
						++possibleOverlapCnt ;
					if ( j - i > longestHits )
						longestHits = j - i  ;
				}
				i = j ;
			}
			// filter based on the repeatability of overlaps.
			if ( possibleOverlapCnt > 10000 )
				novelMinHitRequired = longestHits / 2 ;
			else if ( possibleOverlapCnt > 1000 )
				novelMinHitRequired = longestHits / 3 ;
			else if ( possibleOverlapCnt > 100 )
				novelMinHitRequired = longestHits / 4 ;
		}

		//if ( novelMinHitRequired > 3 )
		//	printf( "novelMinHitRequired=%d\n", novelMinHitRequired ) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].strand != hits[i].strand || hits[j].indexHit.idx != hits[i].indexHit.idx )
					break ;

			int minHitRequired = novelMinHitRequired ;
			if ( seqs[ hits[i].indexHit.idx].isRef )
				minHitRequired = refMinHitRequired ;
		
			/*if ( filter == 1 && readLen > 0 )
			{
				int readStart = hits[i].readOffset, readEnd = hits[j - 1].readOffset + kmerLength - 1 ;
				int seqStart = seqs[ hits[j - 1].indexHit.idx ].consensusLen, seqEnd = -1 ;

				for ( k = i ; k < j ; ++k )
				{
					if ( hits[k].indexHit.offset < seqStart )
						seqStart = hits[k].indexHit.offset ;
					if ( hits[k].indexHit.offset > seqEnd )
						seqEnd = hits[k].indexHit.offset ;
				}
				seqEnd += kmerLength - 1 ;

				int leftOverhangSize = MIN( readStart, seqStart ) ;
				int rightOverhangSize = MIN( readLen - 1 - readEnd, 
						seqs[ hits[i].indexHit.idx ].consensusLen - 1 - seqEnd ) ;
			
				if ( leftOverhangSize > 2 * hitLenRequired || rightOverhangSize > 2 * hitLenRequired )
				{
					i = j ;
					continue ;
				}
			}*/

			//[i,j) holds the hits onto the same seq on the same strand.	
			if ( j - i < minHitRequired )
			{
				i = j ;
				continue ;
			}

			hitCoordDiff.Clear() ;
			for ( k = i ; k < j ; ++k )
			{
				struct _pair nh ;
				nh.a = k ;
				nh.b = hits[k].readOffset - hits[k].indexHit.offset ;
				hitCoordDiff.PushBack( nh ) ;
			}
			std::sort( hitCoordDiff.BeginAddress(), hitCoordDiff.EndAddress(), CompSortPairBInc ) ;

			// Pick the best concordant hits.
			int s, e ;
			int adjustRadius = radius ;
			if ( !seqs[ hits[i].indexHit.idx ].isRef )
				adjustRadius = 0 ;

			for ( s = 0 ; s < j - i ; )
			{
				int diffSum = 0 ;
				for ( e = s + 1 ; e < j - i ; ++e )
				{
					int diff = hitCoordDiff[e].b - hitCoordDiff[e - 1].b ;
					if ( diff < 0 )
						diff = -diff ;
					
					if ( diff > adjustRadius ) 
						break ;

					diffSum += diff ; 
				}
				//printf( "%d %d: %d %d\n", i, j, s, e ) ;


				if ( e - s < minHitRequired 
					|| ( e - s ) * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}

				// [s, e) holds the candidate in the array of hitCoordDiff 
				concordantHitCoord.Clear() ;
				for ( k = s ; k < e ; ++k )
				{
					struct _pair nh ;
					int hitId = hitCoordDiff[k].a ; 
					nh.a = hits[ hitId ].readOffset ;
					nh.b = hits[ hitId ].indexHit.offset ;
					concordantHitCoord.PushBack( nh ) ;
				}
				std::sort( concordantHitCoord.BeginAddress(), concordantHitCoord.EndAddress(), CompSortPairBInc ) ;
				//for ( k = 0 ; k < e - s ; ++k )	
				//	printf( "%d (%d-%d): %d %s %d %d\n", i, s, e, hits[i].indexHit.idx, seqs[ hits[i].indexHit.idx ].name, concordantHitCoord[k].a, concordantHitCoord[k].b ) ;

				// Compute the longest increasing subsequence.
				hitCoordLIS.Clear() ;
				int lisSize = LongestIncreasingSubsequence( concordantHitCoord, hitCoordLIS ) ; 
				if ( lisSize * kmerLength < hitLenRequired )
				{
					s = e ;
					continue ;
				}

				// Rebuild the hits.
				finalHits.Clear() ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _hit nh = hits[i];
					nh.readOffset = hitCoordLIS[k].a ;
					nh.indexHit.offset = hitCoordLIS[k].b ;
					//printf( "%d: %d %d %d %d\n", i, nh.readOffset, nh.indexHit.idx, nh.indexHit.offset, nh.strand ) ;
					finalHits.PushBack( nh ) ;
				}

				int hitLen = GetTotalHitLengthOnRead ( finalHits ) ;
				if ( hitLen < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				else if ( GetTotalHitLengthOnSeq( finalHits ) < hitLenRequired )
				{
					s = e ;
					continue ;
				}
				
				struct _overlap no ;
				no.seqIdx = hits[i].indexHit.idx ;
				no.readStart = finalHits[0].readOffset ;
				no.readEnd = finalHits[ lisSize - 1 ].readOffset + kmerLength - 1 ;
				no.strand = finalHits[0].strand ;
				no.seqStart = finalHits[0].indexHit.offset ;
				no.seqEnd = finalHits[ lisSize - 1 ].indexHit.offset + kmerLength - 1 ;
				no.matchCnt = 2 * hitLen ;
				no.similarity = 0 ;

				if ( !seqs[ no.seqIdx ].isRef && hitLen * 2 < no.seqEnd - no.seqStart + 1 )
				{
					s = e ; 
					continue ;
				}
				no.hitCoords = new SimpleVector<struct _pair> ;
				no.hitCoords->Reserve( lisSize ) ;
				for ( k = 0 ; k < lisSize ; ++k )
				{
					struct _pair nh ;
					nh.a = finalHits[k].readOffset ;
					nh.b = finalHits[k].indexHit.offset ;
					no.hitCoords->PushBack( nh ) ;
				}
				overlaps.push_back( no ) ;

				s = e ;
			} // iterate through concordant hits.
			i = j ;
		}
		return overlaps.size() ;
	}
	
	// Find the overlaps from hits if it possibly span the CDR3 region and anchor paritally on V and J gene 
	int GetVJOverlapsFromHits( SimpleVector<struct _hit> &hits, std::vector<struct _overlap> &overlaps )
	{
		int i, j, k ;
		SimpleVector<struct _hit> VJhits ; 		
		
		int hitSize = hits.Size() ;
		
		// Filter hits that are out of VJ junction region.
		VJhits.Reserve( hitSize ) ;
		for ( i = 0 ; i < hitSize ; ++i )
		{
			int seqIdx = hits[i].indexHit.idx ;
			if ( !seqs[ seqIdx ].isRef )
				continue ;

			if ( seqs[ seqIdx ].name[3] == 'V' && hits[i].indexHit.offset >= seqs[ seqIdx ].consensusLen - 31 )
			{
				VJhits.PushBack( hits[i] ) ;
			}
			else if ( seqs[ seqIdx ].name[3] == 'J' && hits[i].indexHit.offset < 31 )
			{
				VJhits.PushBack( hits[i] ) ;
			}
		}

		GetOverlapsFromHits( VJhits, 17, 0, overlaps ) ;
		
		// Extract the best VJ pair 
		int overlapCnt = overlaps.size() ;
		int maxMatchCnt = 0 ;
		int tagi = 0, tagj = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			for ( j = i + 1 ; j < overlapCnt ; ++j )
			{
				int seqIdxI = overlaps[i].seqIdx ;
				int seqIdxJ = overlaps[j].seqIdx ;

				if ( seqs[ seqIdxI ].name[0] != seqs[ seqIdxJ ].name[0] ||
					seqs[ seqIdxI ].name[1] != seqs[ seqIdxJ ].name[1] ||
					seqs[ seqIdxI ].name[2] != seqs[ seqIdxJ ].name[2] ||
					seqs[ seqIdxI ].name[3] == seqs[ seqIdxJ ].name[3] )
					continue ;			

				if ( seqs[ seqIdxI ].name[3] == 'V' )
				{
					if ( overlaps[i].readStart > overlaps[j].readStart )
						continue ;
				}
				else 
				{
					if ( overlaps[i].readStart < overlaps[j].readStart )
						continue ;
				}

				if ( overlaps[i].matchCnt + overlaps[j].matchCnt > maxMatchCnt )
				{
					maxMatchCnt = overlaps[i].matchCnt + overlaps[j].matchCnt ;
					tagi = i ;
					tagj = j ;
				}
			}
		}

		if ( maxMatchCnt == 0 )
		{
			int size = overlaps.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				overlaps[i].hitCoords->Release() ;
				delete overlaps[i].hitCoords ;
				overlaps[i].hitCoords = NULL ;
			}
			overlaps.clear() ;
			return 0 ;
		}
		int size = overlaps.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( i == tagi || i == tagj )
				continue ;

			overlaps[i].hitCoords->Release() ;
			delete overlaps[i].hitCoords ;
			overlaps[i].hitCoords = NULL ;
		}


		std::vector<struct _overlap> ret ;
		ret.push_back( overlaps[ tagi ] ) ;
		ret.push_back( overlaps[ tagj ] ) ;
		
		overlaps = ret ;

		return 2 ;
	}

	// Extend the overlap to include the overhang parts and filter the overlaps if the overhang does not match well.
	// return: whether this is a valid extension or not
	int ExtendOverlap( char *r, int len, struct _seqWrapper &seq, char *align, struct _overlap &overlap, struct _overlap &extendedOverlap )
	{
		// Check whether the overhang part is compatible with each other or not.
		// Extension to 5'-end ( left end )
		int matchCnt, mismatchCnt, indelCnt ;
		int leftOverhangSize = MIN( overlap.readStart, overlap.seqStart ) ;
		int ret = 1 ;

		//AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqStart - leftOverhangSize, 
		AlignAlgo::GlobalAlignment_PosWeight( seq.posWeight.BeginAddress() + overlap.seqStart - leftOverhangSize, 
				leftOverhangSize, 
				r + overlap.readStart - leftOverhangSize, leftOverhangSize, align ) ;

		GetAlignStats( align, false, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
		{
			leftOverhangSize = 0 ;
			ret = 0 ;
		}
		// Extension to 3'-end ( right end )
		int rightOverhangSize = MIN( len - 1 - overlap.readEnd, seq.consensusLen - 1 - overlap.seqEnd ) ;

		//AlignAlgo::GlobalAlignment( seq.consensus + overlap.seqEnd + 1, 
		AlignAlgo::GlobalAlignment_PosWeight( seq.posWeight.BeginAddress() + overlap.seqEnd + 1, 
				rightOverhangSize,
				r + overlap.readEnd + 1, rightOverhangSize, align ) ;
		GetAlignStats( align, true, matchCnt, mismatchCnt, indelCnt ) ;
		if ( indelCnt > 0 )
		{
			rightOverhangSize = 0 ;
			ret = 0 ;
		}

		int mismatchThreshold = 2 ;
		if ( leftOverhangSize >= 2 )
			++mismatchThreshold ;
		if ( rightOverhangSize >= 2 )
			++mismatchThreshold ;
		
		if ( mismatchCnt > mismatchThreshold && (double)mismatchCnt / ( leftOverhangSize + rightOverhangSize ) > 1.5 / kmerLength ) 
			ret = 0 ;

		extendedOverlap.seqIdx = overlap.seqIdx ;
		extendedOverlap.readStart = overlap.readStart - leftOverhangSize ;
		extendedOverlap.readEnd = overlap.readEnd + rightOverhangSize ;
		extendedOverlap.seqStart = overlap.seqStart - leftOverhangSize ;
		extendedOverlap.seqEnd = overlap.seqEnd + rightOverhangSize ;
		extendedOverlap.strand = overlap.strand ;	
		extendedOverlap.matchCnt = 2 * matchCnt + overlap.matchCnt ;
		extendedOverlap.similarity = (double)( 2 * matchCnt + overlap.matchCnt ) / 
			( extendedOverlap.readEnd - extendedOverlap.readStart + 1 + extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ) ;	
		
		if ( ( seqs[ extendedOverlap.seqIdx ].isRef && extendedOverlap.similarity < refSeqSimilarity ) 
			|| ( !seqs[ extendedOverlap.seqIdx ].isRef && extendedOverlap.similarity < novelSeqSimilarity ) )
		{
			extendedOverlap = overlap ;
			ret = 0 ;
		}

		/*if ( !seqs[ extendedOverlap.seqIdx ].isRef && 
			extendedOverlap.readEnd - extendedOverlap.readStart + 1 - extendedOverlap.matchCnt / 2 >= 5 )
		{
			// Exceed the number of mismatch allowed
			extendedOverlap = overlap ;
			ret = 0 ;
		}*/
		return ret ;
	}

	// Obtain the overlaps, each overlap further contains the hits induce the overlap. 
	// readType: 0(default): sqeuencing read. 1:seqs, no need filter.
	// Return: the number of overlaps.
	int GetOverlapsFromRead( char *read, int readType, std::vector<struct _overlap> &overlaps, 
		SimpleVector<bool> *puse = NULL )
	{
		int i, j, k ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return -1 ;

		SimpleVector<struct _hit> hits ;		    	

		KmerCode kmerCode( kmerLength ) ;
		KmerCode prevKmerCode( kmerLength ) ;

		// Locate the hits from the same-strand case.
		//int skipLimit = 3 ;
		int skipLimit = kmerLength / 2 ; 
		
		int skipCnt = 0 ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;
		
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
			{
				SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

				int size = indexHit.Size() ;
				if ( size >= 100 )
				{
					if ( skipCnt < skipLimit )
					{
						++skipCnt ;
						continue ;
					}
				}
				
				skipCnt = 0 ;

				for ( j = 0 ; j < size ; ++j )
				{
					struct _hit nh ;
					nh.indexHit = indexHit[j] ;
					nh.readOffset = i - kmerLength + 1 ;
					nh.strand = 1 ;
					if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
						continue ;
					hits.PushBack( nh ) ;
				}
			}

			prevKmerCode = kmerCode ;
		}
		// Locate the hits from the opposite-strand case.
		char *rcRead =  new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( rcRead[i] ) ;
		
		skipCnt = 0 ; 
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( rcRead[i] ) ;
			if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
			{
				SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

				int size = indexHit.Size() ;

				if ( size >= 100 )
				{
					if ( skipCnt < skipLimit )
					{
						++skipCnt ;
						continue ;
					}
				}
				skipCnt = 0 ;

				for ( j = 0 ; j < size ; ++j )
				{
					struct _hit nh ;
					nh.indexHit = indexHit[j] ;
					nh.readOffset = i - kmerLength + 1 ;
					nh.strand = -1 ;
					if ( puse != NULL && !puse->Get( indexHit[j].idx ) )
						continue ;
					
					hits.PushBack( nh ) ;
					if ( seqs[indexHit[j].idx].name == NULL)
					{
						printf( "%d %d\n", indexHit[j].idx, indexHit[j].offset ) ;
					}
					assert( seqs[indexHit[j].idx].name != NULL ) ;
				}
			}

			prevKmerCode = kmerCode ;
		}
		delete[] rcRead ;

		//printf( "hitsize=%d; %d\n", hits.Size(), kmerLength ) ;
		// Find the overlaps.
		// Sort the hits
		if ( hits.Size() > 2 * seqs.size() ) 
		{
			// Bucket sort.
			int hitCnt = hits.Size() ;
			int seqCnt = seqs.size() ;
			SimpleVector<struct _hit> *buckets[2] ;
			buckets[0] = new SimpleVector<struct _hit>[seqCnt] ;
			buckets[1] = new SimpleVector<struct _hit>[seqCnt] ;

			for ( i = 0 ; i < hitCnt ; ++i )
			{
				int tag = hits[i].strand == 1 ? 1 : 0 ;
				buckets[tag][ hits[i].indexHit.idx ].PushBack( hits[i] ) ;
			}
			
			hits.Clear() ;
			for ( k = 0 ; k <= 1 ; ++k )
			{
				for ( i = 0 ; i < seqCnt ; ++i )
				{
					hits.PushBack( buckets[k][i] ) ;
				}
			}

			delete[] buckets[0] ;
			delete[] buckets[1] ;
		}
		else
			std::sort( hits.BeginAddress(), hits.EndAddress() ) ;
		//for ( struct _hit *it = hits.BeginAddress() ; it != hits.EndAddress() ; ++it )
		//	printf( "- %d %d %d %d\n", it->readOffset, it->indexHit.idx, it->indexHit.offset, it->strand ) ;

		//int hitLenRequired = 31 ;
		int filterHits = 0 ;
		if ( readType == 0 )
		{
			//hitLenRequired = ( len / 3 < hitLenRequired ? hitLenRequired : ( len / 3 ) ) ;
			filterHits = 1 ;
		}
		int overlapCnt = GetOverlapsFromHits( hits, hitLenRequired, filterHits, overlaps ) ;

		//for ( i = 0 ; i < overlapCnt ; ++i )
		//	printf( "%d: %d %s %d. %d %d %d %d\n", i, overlaps[i].seqIdx,seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd ) ; 
		// Determine whether we want to add this reads by looking at the quality of overlap
		if ( overlapCnt == 0 )
		{
			overlapCnt = GetVJOverlapsFromHits( hits, overlaps ) ;
			if ( overlapCnt == 0 )
				return 0 ;
		}
		// Filter out overlaps that is not a real overlap. 
		/*k = 0 ; 
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			int seqLen = strlen( seqs[ overlaps[i].seqIdx ].consensus ) ; 
			// TODO: adjust the boundary effect for V, J gene.
			if ( ( overlaps[i].readStart - radius >= 0 && overlaps[i].readEnd + radius > len - 1 
						&& overlaps[i].seqStart - radius < 0 ) // the last part of the read overlaps with the index
					|| ( overlaps[i].readStart - radius < 0 && overlaps[i].readEnd + radius <= len - 1 
						&& overlaps[i].seqEnd + radius > seqLen - 1 ) // the first part of the read overlaps with the index
					|| ( overlaps[i].readStart - radius < 0 && overlaps[i].readEnd + radius > len - 1 ) // the read is contained.
			   )
			{
				overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;*/

		// Since the seqs are all from the same strand, the overlaps should be on the same strand.
		std::sort( overlaps.begin(), overlaps.end() ) ;
		
		if ( readType == 0 )
		{
			k = 1 ;
			for ( i = 1 ; i < overlapCnt ; ++i )
			{
				if ( overlaps[i].strand != overlaps[0].strand )
				{
					delete overlaps[i].hitCoords ;
					overlaps[i].hitCoords = NULL ;
					continue ;		
				}
				if ( i != k )
					overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		else
		{
			k = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				if ( overlaps[i].strand != 1 )
				{
					delete overlaps[i].hitCoords ;
					overlaps[i].hitCoords = NULL ;
					continue ;		
				}
				if ( i != k )
					overlaps[k] = overlaps[i] ;
				++k ;
			}
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		/*for ( i = 1 ; i < overlapCnt ; ++i )
		{
			if ( overlaps[i].strand != overlaps[i - 1].strand )
			{
				overlaps.clear() ;
				for ( i = 0 ; i < overlapCnt ; ++i )
					delete overlaps[i].hitCoords ;
				return 0 ;
			}
		}*/

		rcRead = new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		

		// Compute similarity overlaps
		//if ( overlaps.size() > 1000 )
		/*int bestNovelOverlap = -1 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef )
				continue ;

			if ( bestNovelOverlap == -1 || overlaps[i] < overlaps[ bestNovelOverlap ] )
				bestNovelOverlap = i ;
		}
		if ( bestNovelOverlap != -1 )
		{
			struct _overlap tmp ;
			tmp = overlaps[0] ;
			overlaps[0] = overlaps[ bestNovelOverlap ] ;
			overlaps[ bestNovelOverlap ] = tmp ;
		}*/
		
		int firstRef = -1 ;
		int bestNovelOverlap = -1 ;
		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *r ;
			if ( overlaps[i].strand == 1 )
				r = read ;
			else
				r = rcRead ;

			SimpleVector<struct _pair> &hitCoords = *overlaps[i].hitCoords ; 	
			int hitCnt = hitCoords.Size() ;
			int matchCnt = 0, mismatchCnt = 0, indelCnt = 0  ;
			double similarity = 1 ;
			
			// Some fast pre-filters
			if ( seqs[ overlaps[i].seqIdx ].isRef )
			{
				if ( firstRef == -1 )
					firstRef = i ;
				/*else
				{
					if ( overlaps[i].matchCnt < 0.9 * overlaps[ firstRef ].matchCnt )
					{
						overlaps[i].similarity = 0 ;  // No need to look into this.
						continue ;
					}
				}*/
			}
			else if ( bestNovelOverlap != -1 && readType == 0 && overlapCnt > 50 )
			{
				// If we already found a perfect match.
				if ( overlaps[ bestNovelOverlap ].readStart == 0 && overlaps[ bestNovelOverlap ].readEnd == len - 1 ) 
				{
					//printf( "fast perfect filter %d: %lf %d %d\n", overlaps[ bestNovelOverlap ].seqIdx,
					//	overlaps[ bestNovelOverlap ].similarity, overlaps[ bestNovelOverlap ].matchCnt, 
					//	overlaps[i].matchCnt ) ;
					if ( overlaps[ bestNovelOverlap ].similarity == 1 )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
					else if ( overlaps[ bestNovelOverlap ].similarity > repeatSimilarity &&
						overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
				}
				
				// Almost fully cover case
				if ( overlaps[ bestNovelOverlap ].readStart + len - 1 - overlaps[ bestNovelOverlap].readEnd < radius )
				{
					if ( overlaps[ bestNovelOverlap].similarity == 1 
						&& overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
					else if ( overlaps[ bestNovelOverlap ].similarity > repeatSimilarity 
						&& overlaps[i].matchCnt < 0.8 * overlaps[ bestNovelOverlap ].matchCnt )
					{
						overlaps[i].similarity = 0 ;
						continue ;
					}
				}
				
				// Directly filter the bad overlaps if it is inside of a novel seq and worse than the best one. 
				if ( overlaps[i].seqStart - overlaps[i].readStart >= radius 
						&& overlaps[i].seqEnd + ( len - 1 - overlaps[i].readEnd ) + radius < 
							seqs[ overlaps[i].seqIdx ].consensusLen 
						&& overlaps[ bestNovelOverlap ].matchCnt > 0.97 * ( 2 * len ) 
						&& overlaps[ bestNovelOverlap ].similarity > repeatSimilarity 
						&& overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					//printf( "fast filter\n" ) ;
					overlaps[i].similarity = 0 ;	
					continue ;
				}

				if ( overlaps[i].matchCnt < 0.4 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}

				if ( overlapCnt > 1000 && overlaps[i].matchCnt < 0.9 * overlaps[ bestNovelOverlap ].matchCnt )
				{
					overlaps[i].similarity = 0 ;
					continue ;
				}

				//printf( "%d %d: %d %d\n", overlapCnt, overlaps[i].seqIdx, overlaps[i].matchCnt, overlaps[ bestNovelOverlap ].matchCnt ) ;
			}


			matchCnt += 2 * kmerLength ;
			char *align = new char[ overlaps[i].readEnd - overlaps[i].readStart + 1 + 
				overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 1] ;
			for ( j = 1 ; j < hitCnt ; ++j )
			{
				if ( hitCoords[j - 1].b - hitCoords[j - 1].a == hitCoords[j].b - hitCoords[j].a )
				{
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ; 
						
						if ( seqs[ overlaps[i].seqIdx ].isRef  )
						{
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;
						}
						else
						{
						//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) ,
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength), 
								align ) ;
						}

						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;

						if ( ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef ) && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}
					}
				}
				else
				{
					if ( radius == 0 || !seqs[ overlaps[i].seqIdx ].isRef )
					{
						similarity = 0 ;
						break ;
					}

					//printf( "%d %d=>%d %d\n", hitCoords[j - 1].a, hitCoords[j - 1].b, hitCoords[j].a, hitCoords[j].b ) ;
					if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 < hitCoords[j].b )
					{
						matchCnt += 2 * ( hitCoords[j].a - hitCoords[j - 1].a ) ; //+ kmerLength ;
						// Make the two kmer hit match on coordinate.
						indelCnt += ( hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ) + 
							( hitCoords[j].a + kmerLength - hitCoords[j - 1].a )  ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 < hitCoords[j].a && 
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += kmerLength + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ( hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) +
							( hitCoords[j].b + kmerLength - hitCoords[j - 1].b ) ) ;
					}
					else if ( hitCoords[j - 1].a + kmerLength - 1 >= hitCoords[j].a &&
						hitCoords[j - 1].b + kmerLength - 1 >= hitCoords[j].b )
					{
						//matchCnt += ( hitCoords[j].a - hitCoords[j - 1].a ) + ( hitCoords[j].b - hitCoords[j - 1].b ) ;
						matchCnt += 2 * MIN( hitCoords[j].a - hitCoords[j - 1].a, hitCoords[j].b - hitCoords[j - 1].b ) ;
						indelCnt += ABS( ( hitCoords[j].a - hitCoords[j].b ) - 
							( hitCoords[j - 1].a - hitCoords[j - 1].b ) ) ;
					}
					else
					{
						matchCnt += 2 * kmerLength ;
						 
						if ( seqs[ overlaps[i].seqIdx ].isRef )
						{
							AlignAlgo::GlobalAlignment( 
								seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						else
						{
							//AlignAlgo::GlobalAlignment( seqs[ overlaps[i].seqIdx ].consensus + hitCoords[j - 1].b + kmerLength,
							AlignAlgo::GlobalAlignment_PosWeight( 
								seqs[ overlaps[i].seqIdx ].posWeight.BeginAddress() + hitCoords[j - 1].b + kmerLength,
								hitCoords[j].b - ( hitCoords[j - 1].b + kmerLength ),
								r + hitCoords[j - 1].a + kmerLength, 
								hitCoords[j].a - ( hitCoords[j - 1].a + kmerLength ) , 
								align ) ;	
						}
						
						int count[3] ;
						GetAlignStats( align, false, count[0], count[1], count[2] ) ;
						matchCnt += 2 * count[0] ;
						mismatchCnt += count[1] ;
						indelCnt += count[2] ;
						if ( !seqs[ overlaps[i].seqIdx ].isRef && indelCnt > 0 )
						{
							similarity = 0 ;
							break ;
						}

					}
				}
			} // for j
			delete[] align ;
			
			//printf( "%d %d %d %lf\n", matchCnt, overlaps[i].seqEnd - overlaps[i].seqStart + 1, overlaps[i].readEnd - overlaps[i].readStart + 1, similarity ) ;
			overlaps[i].matchCnt = matchCnt ;
			if ( similarity == 1 )
				overlaps[i].similarity = (double)matchCnt / ( overlaps[i].seqEnd - overlaps[i].seqStart + 1 + 
								overlaps[i].readEnd - overlaps[i].readStart + 1 ) ;
			else
				overlaps[i].similarity = 0 ;
			
			if ( IsOverlapLowComplex( r, overlaps[i]) )
				overlaps[i].similarity = 0 ;
			

			if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity > 0 )
			{
				if ( bestNovelOverlap == -1 || overlaps[i] < overlaps[ bestNovelOverlap ] ) // the less than means has higher priority
				{
					bestNovelOverlap = i ;
				}
			}
			/*if ( overlaps[i].similarity > 0 )
			{
				printf( "%d: %d %d %d %d %d %lf\n", matchCnt, overlaps[i].seqIdx, overlaps[i].readStart, overlaps[i].readEnd, 
							overlaps[i].seqStart, overlaps[i].seqEnd, similarity ) ;
			}
			assert( overlaps[i].similarity <= 1 ) ;*/
		} // for i
		delete[] rcRead ;

		// Release the memory for hitCoords.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			overlaps[i].hitCoords->Release() ;
			delete overlaps[i].hitCoords ;
			overlaps[i].hitCoords = NULL ;
		}

		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < refSeqSimilarity )
				continue ;
			else if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].similarity < novelSeqSimilarity )
				continue ;

			//printf( "%lf\n", overlaps[i].similarity ) ;
			overlaps[k] = overlaps[i] ;
			++k ;
		}
		overlaps.resize( k ) ;
		overlapCnt = k ;
		

		return overlapCnt ;
	}
	
	// Figure out whether a seq is a (almost) substring of another seq.
	int BuildSeqSubstringRelation( std::vector<struct _overlap> &subsetOf )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		for ( k = 0 ; k < seqCnt ; ++k )
			subsetOf[k].seqIdx = -1 ;
		
		std::map<int, int> seqHitCnt ;
		SimpleVector<struct _pair> firstSeqHit ; 
		firstSeqHit.ExpandTo( seqCnt ) ;

		for ( k = 0 ; k < seqCnt ; ++k )
		{
			if ( seqs[k].consensus == NULL )
				continue ;

			char *consensus = seqs[k].consensus ;
			int len = seqs[k].consensusLen ;
			if ( len < kmerLength )
				return -1 ;

			//SimpleVector<struct _hit> hits ;		    	

			KmerCode kmerCode( kmerLength ) ;
			KmerCode prevKmerCode( kmerLength ) ;

			//int skipLimit = 3 ;
			int skipLimit = kmerLength / 2 ; 

			int skipCnt = 0 ;
			int hitCnt ;
			for ( i = 0 ; i < kmerLength - 1 ; ++i )
				kmerCode.Append( consensus[i] ) ;
			seqHitCnt.clear() ;
			hitCnt = 0 ;
			for ( ; i < len ; ++i )
			{
				kmerCode.Append( consensus[i] ) ;
				if ( i == kmerLength || !prevKmerCode.IsEqual( kmerCode ) )
				{
					SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

					int size = indexHit.Size() ;
					if ( size >= 100 )
					{
						if ( skipCnt < skipLimit )
						{
							++skipCnt ;
							continue ;
						}
					}

					skipCnt = 0 ;

					for ( j = 0 ; j < size ; ++j )
					{
						struct _hit nh ;
						if ( indexHit[j].idx == k )
							continue ;
												
						if ( seqHitCnt.find( indexHit[j].idx ) != seqHitCnt.end() )
						{
							if ( hitCnt >= 50 && seqHitCnt[ indexHit[j].idx ] < hitCnt * 0.5 )
							{
								seqHitCnt.erase( indexHit[j].idx ) ;
							}
							else
								++seqHitCnt[ indexHit[j].idx ] ;
						}
						else if ( hitCnt < 50 )
						{
							seqHitCnt[ indexHit[j].idx ] = 1 ;
							firstSeqHit[ indexHit[j].idx ].a = i - kmerLength + 1 ;
							firstSeqHit[ indexHit[j].idx ].b = indexHit[j].offset ;
						}
					}
					++hitCnt ;
				}

				prevKmerCode = kmerCode ;
			}
			
			for ( std::map<int, int>::iterator it = seqHitCnt.begin() ; it != seqHitCnt.end() ; ++it )
			{
				if ( it->second < hitCnt * 0.75 )
					continue ;
				int seqIdx = it->first ; 
				// Test whether k is a substring of seqIdx
				if ( firstSeqHit[ seqIdx ].b - firstSeqHit[ seqIdx ].a < 0 )
					continue ;

				int start = firstSeqHit[ seqIdx ].b - firstSeqHit[ seqIdx ].a ;
				if ( start + seqs[k].consensusLen - 1 >= seqs[seqIdx].consensusLen )
					continue ;
				int matchCnt = 0 ;
				int mismatchCnt = 0 ;
				int l ;
				for ( j = 0, l = start ; j < seqs[k].consensusLen ; ++j, ++l )
				{
					if ( seqs[k].consensus[j] != seqs[seqIdx].consensus[l] )
						++mismatchCnt ;
					else
						++matchCnt ;

					if ( mismatchCnt >= 2 )
						break ;
				}
				//char *p = strstr( seqs[ it->first ].consensus, consensus ) ;
				//printf( "test %d\n%s\n%s\n", mismatchCnt, seqs[k].consensus, seqs[ seqIdx ].consensus  ) ;
				if ( mismatchCnt < 2 ) // some mismatch are allowed because we allow mismatch in the overlaps of branch graph.
				{
					subsetOf[k].seqIdx = it->first ;
					subsetOf[k].readStart = 0 ;
					subsetOf[k].readEnd = seqs[k].consensusLen - 1 ;
					subsetOf[k].seqStart = start ;
					subsetOf[k].seqEnd = subsetOf[k].seqStart + seqs[k].consensusLen - 1 ; 
					break ;
				}
			}
		}

		return 0 ;
	}


	// adj is the adjacent list for each seq. The size of the adj array is seqCnt.
	int BuildSeqOverlapGraph( int overlapLength, std::vector<struct _overlap> *adj )
	{
		// Build overlap graph.
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int maxLen = 0 ;
		char *align ;
		for ( i = 0 ; i < seqCnt ; ++i )
			if ( seqs[i].consensusLen > maxLen )
				maxLen = seqs[i].consensusLen ;
		align = new char[2 * maxLen + 2] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL ) 
				continue ;

			std::vector<struct _overlap> overlaps ;
			int overlapCnt ;

			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, 1, overlaps ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx  )
					continue ;
				struct _overlap extendedOverlap ;

				if ( ExtendOverlap( seqs[i].consensus, seqs[i].consensusLen, seqs[ overlaps[j].seqIdx ], 
					align, overlaps[j], extendedOverlap ) == 1 ) 
				{
					if ( extendedOverlap.readEnd - extendedOverlap.readStart + 1 >= overlapLength && 
						( extendedOverlap.readStart > 0 || // i before j
						( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == seqs[i].consensusLen - 1 )// i contained in j 
						) ) 
					{
						if ( ( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == seqs[i].consensusLen - 1 ) 
							&& ( extendedOverlap.seqStart == 0 && 
								extendedOverlap.seqEnd == seqs[ overlaps[j].seqIdx ].consensusLen - 1 ) )
						{
							// The case of two almost identical seqs.
							if ( i < overlaps[j].seqIdx )
								continue ;
						}
						adj[i].push_back( extendedOverlap ) ;	
					}
				}
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
			std::sort( adj[i].begin(), adj[i].end() ) ;
		delete[] align ;
		return 1 ;
	}
	
	// Only build the graph for the seqs with id marked true in "use" array.
	int BuildBranchGraph( std::vector<struct _overlap> *adj, int leastOverlapLen, SimpleVector<bool> &use )
	{
		// Build overlap graph.
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int maxLen = 0 ;
		char *align ;
		for ( i = 0 ; i < seqCnt ; ++i )
			if ( seqs[i].consensusLen > maxLen )
				maxLen = seqs[i].consensusLen ;
		align = new char[2 * maxLen + 2] ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL ) 
				continue ;
			if ( use[i] == false )
				continue ;

			std::vector<struct _overlap> overlaps ;
			int overlapCnt ;

			overlapCnt = GetOverlapsFromRead( seqs[i].consensus, 1, overlaps, &use ) ;
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( overlaps[j].strand == -1 ) // Note that all the seqs are from 5'->3' on its strand.
					continue ;

				if ( i == overlaps[j].seqIdx || use[ overlaps[j].seqIdx ] == false )
					continue ;
				struct _overlap extendedOverlap ;
				int addMatchCnt = 0 ;				
				// Locally extend the overlap. In this case, the overlap actually just means share.
				// Allow "radius" overhang.
				int seqIdx = overlaps[j].seqIdx ;
				int a, b, k ;
				int matchCnt = 0 ;
				int rightExtend = 0;
				int rightExtendMatchCnt = 0 ; 
				for ( k = 1, a = overlaps[j].readEnd + 1, b = overlaps[j].seqEnd + 1 ; 
					a < seqs[i].consensusLen && b < seqs[ seqIdx ].consensusLen ; ++a, ++b, ++k )
				{
					if ( IsPosWeightCompatible( seqs[i].posWeight[a], seqs[ seqIdx ].posWeight[b] ) )
						++matchCnt ;

					if ( matchCnt > k * 0.75 )
					{
						rightExtendMatchCnt = 2 * matchCnt ;
						rightExtend = k ;
					}
				}

				matchCnt = 0 ;
				int leftExtend = 0 ;
				int leftExtendMatchCnt = 0 ;
				for ( k = 1, a = overlaps[j].readStart - 1, b = overlaps[j].seqStart - 1 ; 
					a >= 0 && b >= 0 ; --a, --b, ++k )
				{
					/*if ( b >= seqs[ seqIdx ].consensusLen )
						fprintf( stderr, "%d %d %d %d\n%s\n%s\n", overlaps[j].readStart, overlaps[j].readEnd,
								overlaps[j].seqStart, overlaps[j].seqEnd,
								seqs[i].consensus, seqs[ seqIdx ].consensus ) ;*/
					if ( IsPosWeightCompatible( seqs[i].posWeight[a], seqs[ seqIdx ].posWeight[b] ) )
						++matchCnt ;

					if ( matchCnt > k * 0.75 )
					{
						leftExtendMatchCnt = 2 * matchCnt ;
						leftExtend = k ;
					}
				}
				
				extendedOverlap = overlaps[j] ;
				extendedOverlap.readStart -= leftExtend ;
				extendedOverlap.seqStart -= leftExtend ;
				extendedOverlap.readEnd += rightExtend ;
				extendedOverlap.seqEnd += rightExtend ;
				extendedOverlap.matchCnt += rightExtendMatchCnt + leftExtendMatchCnt ;
				extendedOverlap.similarity = (double)extendedOverlap.matchCnt / 
					( extendedOverlap.readEnd - extendedOverlap.readStart + 1 +
					  extendedOverlap.seqEnd - extendedOverlap.seqStart + 1 ) ;
				
				if ( extendedOverlap.readEnd - extendedOverlap.readStart + 1 < leastOverlapLen )
					continue ;


				if ( extendedOverlap.similarity >= repeatSimilarity )
				{
					adj[i].push_back( extendedOverlap ) ;
#ifdef DEBUG
					printf( "branch %d %d: %d %d %d %d %d %lf\n", i, j, extendedOverlap.seqIdx, 
							extendedOverlap.readStart, extendedOverlap.readEnd,
							extendedOverlap.seqStart, extendedOverlap.seqEnd,
							extendedOverlap.similarity ) ;
#endif				
				}	
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
			std::sort( adj[i].begin(), adj[i].end() ) ;
		delete[] align ;
		return 1 ;
	}


	void UpdatePosWeightFromRead( SimpleVector<struct _posWeight> &posWeight, int offset, char *read )
	{
		int i ;
		for ( i = 0 ; read[i] ; ++i )
		{
			if ( read[i] != 'N' )
				++posWeight[i + offset].count[ nucToNum[ read[i] - 'A' ] ] ;
		}
	}
public:
	SeqSet( int kl ) 
	{
		kmerLength = kl ;
		radius = 10 ;
		hitLenRequired = 31 ;

		novelSeqSimilarity = 0.9 ;
		refSeqSimilarity = 0.75 ; 
		repeatSimilarity = 0.95 ; 

		prevAddInfo.readStart = -1 ;
	}
	~SeqSet() 
	{
		int size ;
		int i ;
		size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].consensus != NULL )
				free( seqs[i].consensus ) ;
			if ( seqs[i].name != NULL )
				free( seqs[i].name ) ;	
		}
	}

	int Size()
	{
		return seqs.size() ;
	}

	int SetRadius( int r )  
	{
		return radius = r ;
	}
	
	void ReverseComplement( char *rcSeq, char *seq, int len )
	{
		int i ;
		for ( i = 0 ; i < len ; ++i )
		{
			if ( seq[len - 1 - i] != 'N' )
				rcSeq[i] = numToNuc[ 3 - nucToNum[seq[len - 1 - i] - 'A'] ];
			else
				rcSeq[i] = 'N' ;
		}
		rcSeq[i] = '\0' ;
	}

	// Input some baseline sequence to match against.
	void InputRefFa( char *filename ) 
	{
		int i, j, k ;
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		KmerCode kmerCode( kmerLength ) ;
		while ( fa.Next() )
		{
			// Insert the kmers 
			struct _seqWrapper ns ;
			ns.name = strdup( fa.id ) ;
			ns.isRef = true ;

			int id = seqs.size() ;
			seqs.push_back( ns ) ;

			struct _seqWrapper &sw = seqs[id] ;
			int seqLen = strlen( fa.seq ) ;
			sw.consensus = strdup( fa.seq ) ;	
			
			
			// Remove "." from IMGT annotation.
			k = 0 ;
			for ( i = 0 ; i < seqLen ; ++i )
				if ( sw.consensus[i] != '.' )
				{
					sw.consensus[k] = sw.consensus[i] ;
					++k ;
				}
			sw.consensus[k] = '\0' ;

			// Use IMGT documented coordinate to infer CDR1,2 coordinate.
			if ( GetGeneType( fa.id ) == 0 && seqLen >= 66 * 3 )
			{	
				// Infer the coordinate for CDR1
				for ( i = 0, k = 0 ; i < 3 * ( 27 - 1 ) ; ++i )
					if ( fa.seq[i] != '.' )		
						++k ;
				sw.info[0].a = k ;
				for ( ; i < 3 * ( 38 ) ; ++i )
					if ( fa.seq[i] != '.' )	
						++k ;
				sw.info[0].b = k - 1 ;

				// Infer the coordinate for CDR2
				for ( ; i < 3 * ( 56 - 1 ) ; ++i )
					if ( fa.seq[i] != '.' )		
						++k ;
				sw.info[1].a = k ;
				for ( ; i < 3 * ( 65 ) ; ++i )
					if ( fa.seq[i] != '.' )	
						++k ;
				sw.info[1].b = k - 1 ;
			}
			else
			{
				sw.info[0].a = sw.info[0].b = -1 ;
				sw.info[1].a = sw.info[1].b = -1 ;
			}

			

			sw.consensusLen = strlen( sw.consensus );	
			seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, sw.consensusLen, id ) ;
		}
	}
	
	void InputNovelFa( char *filename ) 
	{
		ReadFiles fa ;
		fa.AddReadFile( filename, false ) ;
		
		while ( fa.Next() )
			InputNovelRead( fa.id, fa.seq, 1 ) ;
	}

	int InputNovelRead( char *id, char *read, int strand )
	{
		struct _seqWrapper ns ;
		ns.name = strdup( id ) ;
		ns.isRef = false ;

		int seqIdx = seqs.size() ;
		seqs.push_back( ns ) ;

		struct _seqWrapper &sw = seqs[ seqIdx ] ;
		int seqLen = strlen( read ) ;
		sw.consensus = strdup( read ) ;
		if ( strand == -1 )
			ReverseComplement( sw.consensus, read, seqLen ) ;
		sw.consensusLen = seqLen ;	
		int i ;
		sw.posWeight.ExpandTo( sw.consensusLen ) ;
		sw.posWeight.SetZero( 0, sw.consensusLen ) ;
		for ( i = 0 ; i < sw.consensusLen ; ++i )
		{
			if ( sw.consensus[i] != 'N' )
				sw.posWeight[i].count[ nucToNum[ sw.consensus[i] - 'A' ] ] = 1 ;
		}
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.BuildIndexFromRead( kmerCode, sw.consensus, seqLen, seqIdx ) ;
	
		sw.minLeftExtAnchor = sw.minRightExtAnchor = 0 ;
		SetPrevAddInfo( seqIdx, 0, seqLen - 1, 0, seqLen - 1, strand ) ;


#ifdef DEBUG
		printf( "add novel seq: %d\n", seqIdx ) ;
#endif
		return seqIdx ;
	}

	// Input sequence from another SeqSet
	void InputSeqSet( const SeqSet &in, bool inputRef )
	{
		int i, k ;
		int seqCnt = in.seqs.size() ;
		KmerCode kmerCode( kmerLength ) ;
		
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( !inputRef && in.seqs[i].isRef )
				continue ;
			if ( in.seqs[i].consensus == NULL )
				continue ;
			
			struct _seqWrapper ns ;
			ns.name = strdup( in.seqs[i].name ) ;
			ns.isRef = in.seqs[i].isRef ;

			int id = seqs.size() ;
			ns.consensus = strdup( in.seqs[i].consensus ) ;
			ns.consensusLen = in.seqs[i].consensusLen ;
			ns.posWeight = in.seqs[i].posWeight ;
			ns.minLeftExtAnchor = in.seqs[i].minLeftExtAnchor ;
			ns.minRightExtAnchor = in.seqs[i].minRightExtAnchor ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, ns.consensusLen, id ) ;
			seqs.push_back( ns ) ;
		}
	}
	
	// Test whether the read share a kmer hit on the seqs.
	bool HasHitInSet( char *read )
	{
		int i ;
		int len = strlen( read ) ;
		if ( len < kmerLength )
			return false ;

		KmerCode kmerCode( kmerLength ) ;

		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( read[i] ) ;
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( read[i] ) ;
			SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 
			int size = indexHit.Size() ;
			if ( size > 0 )
				return true ;
		}
		
		// Locate the hits from the opposite-strand case.
		char *rcRead =  new char[len + 1] ;
		ReverseComplement( rcRead, read, len ) ;		
		kmerCode.Restart() ;
		for ( i = 0 ; i < kmerLength - 1 ; ++i )
			kmerCode.Append( rcRead[i] ) ;
		
		for ( ; i < len ; ++i )
		{
			kmerCode.Append( rcRead[i] ) ;
			SimpleVector<struct _indexInfo> &indexHit = *seqIndex.Search( kmerCode ) ; 

			int size = indexHit.Size() ;
			if ( size > 0 )
			{
				delete[] rcRead ;
				return true ;
			}
		}

		delete[] rcRead ;
		return false ;
	}

	// Compute the length of hit from the read, take the overlaps of kmer into account 
	int GetTotalHitLengthOnRead( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;
		//for ( i = 0 ; i < hitSize ; ++i )
		//	printf( "%d %d\n", i, hits[i].readOffset) ;
		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
			{
				if ( hits[j].readOffset > hits[j - 1].readOffset + kmerLength - 1 )	
					break ;
			}

			ret += hits[j - 1].readOffset - hits[i].readOffset + kmerLength ;

			i = j ;
		}
		return ret ;
	}

	int GetTotalHitLengthOnSeq( SimpleVector<struct _hit> &hits )
	{
		int hitSize = hits.Size() ;
		int i, j ;
		int ret = 0 ;

		for ( i = 0 ; i < hitSize ; )
		{
			for ( j = i + 1 ; j < hitSize ; ++j )
				if ( hits[j].indexHit.offset > hits[j - 1].indexHit.offset + kmerLength - 1 )
					break ;
			ret += hits[j - 1].indexHit.offset - hits[i].indexHit.offset + kmerLength ;
			i = j ;
		}
		return ret ;
	}

	// b comes after a, test wheter their gene names is compatible
	bool IsNameCompatible( char *a, char *b )
	{
		int i, j ;
		int maxA = -1 ;
		int minB = 10 ;
		for ( i = 0 ; a[i] ; )	
		{
			if ( a[i] == '+' )
			{
				++i ;
				continue ;
			}

			for ( j = i ; a[j] && a[j] != '+' ; ++j )
				;
			char tmp = a[j] ;
			a[j] = '\0' ;
			int gt = GetGeneType( a + i ) ;
			if ( gt > maxA )
				maxA = gt ;
			a[j] = tmp ;

			i = j ;
		}
		
		for ( i = 0 ; b[i] ; )	
		{
			if ( b[i] == '+' )
			{
				++i ;
				continue ;
			}

			for ( j = i ; b[j] && b[j] != '+' ; ++j )
				;
			char tmp = b[j] ;
			b[j] = '\0' ;
			int gt = GetGeneType( b + i ) ;
			if ( gt < minB )
				minB = gt ;
			b[j] = tmp ;

			i = j ;
		}

		if ( maxA <= minB )
			return true ;
		else
			return false ;
	}
	

	// Test whether a read can from the index and update the index.
	// If it is a candidate, but is quite different from the one we stored, we create a new poa for it.
	// Return: the index id in the set. 
	//	   -1: not add. -2: only overlapped with novel seq and could not be extended.
	int AddRead( char *read, char *geneName, double similarityThreshold )
	{
		//printf( "%s\n", read ) ;
		int i, j, k ;
		int len = strlen( read ) ;

		SetPrevAddInfo( -1, -1, -1, -1, -1, 0 ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
		
		overlapCnt = GetOverlapsFromRead( read, 0, overlaps ) ;
				
		if ( overlapCnt <= 0 )
			return -1 ;

#ifdef DEBUG
		printf( "geneName: %s\n", geneName ) ;
#endif
		if ( geneName[0] != '\0' )
		{
			k = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				for ( j = 0 ; j < 3 ; ++j )
					if ( seqs[ overlaps[i].seqIdx ].name[j] != geneName[j] )
						break ;
				if ( j == 3 )
				{
					overlaps[k] = overlaps[i] ;
					++k ;
				}
				/*else
				{
					printf( "haha: %s %s\n%s\n", geneName, seqs[ overlaps[i].seqIdx ].name, read ) ;
				}*/
			}
			overlaps.resize( k ) ;
			overlapCnt = k ;
			
			if ( overlapCnt <= 0 )
				return -1 ;
		}

		std::sort( overlaps.begin(), overlaps.end() ) ;

#ifdef DEBUG
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			printf( "%d: %d %d %s. %d. %d %d %d %d. %lf.\n", i, overlaps[i].seqIdx, seqs[ overlaps[i].seqIdx ].consensusLen, 
					seqs[ overlaps[i].seqIdx ].name, overlaps[i].strand, 
					overlaps[i].readStart, overlaps[i].readEnd, overlaps[i].seqStart, overlaps[i].seqEnd, 
					overlaps[i].similarity ) ; 
			printf( "%s\n", seqs[  overlaps[i].seqIdx ].consensus ) ;
			//if ( !seqs[ overlaps[i].seqIdx ].isRef )
			/*if ( !strcmp( read, "TTACTGTAATATACGATATTTTGACTGGTTATTAAGAGGCGACCCAAGAATCAATACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCT" ) )
			{
				printf( " %s\n",seqs[  overlaps[i].seqIdx ].consensus ) ;
				if ( i == 1 )
				{
					printf( "%d %d\n", seqs[ overlaps[i].seqIdx ].posWeight[103].count[0], 
						seqs[ overlaps[i].seqIdx ].posWeight[103].count[2] ) ;
				}
			}*/
		}
		fflush( stdout ) ;	
#endif		
		// If the read only overlaps with the reference, we will add that to the seq.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( !seqs[ overlaps[i].seqIdx ].isRef )
				break ;
		}
		
		struct _overlap *extendedOverlaps = new struct _overlap[ overlapCnt ];
		struct _overlap *failedExtendedOverlaps = new struct _overlap[ overlapCnt ];
		k = 0 ;
		int ret = -1 ;
		bool addNew = true ;
		int failedExtendedOverlapsCnt = 0 ;

		if ( i < overlapCnt )
		{
			// Incorporate to existing sequence.
			char *align = new char[3 * len] ;
			char *rcRead = strdup( read ) ;
			ReverseComplement( rcRead, read, len ) ;
			
			char *r ;
			int readInConsensusOffset = 0 ;
			int seqIdx ;
			int tag ;
			bool sortExtendedOverlaps = true ;

			if ( overlaps[0].strand == 1 )
				r = read ;
			else
				r = rcRead ;

#ifdef DEBUG
			if ( overlaps[0].strand == -1 )
				printf( "rc: %s\n", r ) ;
#endif
			/*int tmpCnt = 0 ;
			for ( i = 0 ; i < overlapCnt ; ++i )
				if ( !seqs[ overlaps[i].seqIdx ].isRef && overlaps[i].readStart == 0 && overlaps[i].readEnd == 74 && overlaps[i].similarity == 1 )
				{
					++tmpCnt ;
				}
			assert ( tmpCnt <= 1 ) ;*/

			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				for ( j = 0 ; j < k ; ++j )
				{
					// If is contained in extended sequence.
					int leftRadius = radius ;
					int rightRadius = radius ;
					if ( extendedOverlaps[j].seqStart == 0 )
						leftRadius = 0 ;
					if ( extendedOverlaps[j].seqEnd == seqs[ extendedOverlaps[j].seqIdx ].consensusLen - 1 )
						rightRadius = 0 ;
					if ( overlaps[i].readStart >= extendedOverlaps[j].readStart - leftRadius  
						&& overlaps[i].readEnd <= extendedOverlaps[j].readEnd + rightRadius )
						break ;
					
					// Some extended is a subset of this one.
					leftRadius = radius ;
					rightRadius = radius ;
					if ( overlaps[i].seqStart == 0 )
						leftRadius = 0 ;
					if ( overlaps[i].seqEnd == seqs[ overlaps[i].seqIdx ].consensusLen - 1 )
						rightRadius = 0 ;
					if ( extendedOverlaps[j].readStart >= overlaps[i].readStart - leftRadius 
						&& extendedOverlaps[j].readEnd <= overlaps[i].readEnd + rightRadius )
						break ;
				}
				
				struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
				if ( j < k || seq.isRef )
					continue ;
				//if ( !strcmp( read, "TTACTGTAATATACGATATTTTGACTGGTTATTAAGAGGCGACCCAAGAATCAATACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCT" ) && i == 1 )
				//	fprintf( stderr, "hi\n" ) ;
				// Only extend the novel seqs.
				if ( ExtendOverlap( r, len, seq, align, overlaps[i], extendedOverlaps[k] ) == 1 )
				{
					// Double check whether there is subset relationship.
					for ( j = 0 ; j < k ; ++j )
					{
						int leftRadius = radius ;
						int rightRadius = radius ;
						if ( extendedOverlaps[j].seqStart == 0 )
							leftRadius = 0 ;
						if ( extendedOverlaps[j].seqEnd == seqs[ extendedOverlaps[j].seqIdx ].consensusLen - 1 )
							rightRadius = 0 ;
						if ( extendedOverlaps[k].readStart >= extendedOverlaps[j].readStart - leftRadius  
								&& extendedOverlaps[k].readEnd <= extendedOverlaps[j].readEnd + rightRadius )
							break ;

						// Some extended is a subset of this one.
						leftRadius = radius ;
						rightRadius = radius ;
						if ( extendedOverlaps[k].seqStart == 0 )
							leftRadius = 0 ;
						if ( extendedOverlaps[k].seqEnd == seqs[ extendedOverlaps[k].seqIdx ].consensusLen - 1 )
							rightRadius = 0 ;
						if ( extendedOverlaps[j].readStart >= extendedOverlaps[k].readStart - radius 
								&& extendedOverlaps[j].readEnd <= extendedOverlaps[k].readEnd + radius )
							break ;
					}
					if ( j < k )
						continue ;

					// Then check whether the extended porition is a subset of matched portion from other overlaps.
					for ( j = 0 ; j < i ; ++j )
					{
						if ( seqs[ overlaps[j].seqIdx ].isRef )
							continue ;

						if ( extendedOverlaps[k].readStart >= overlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= overlaps[j].readEnd )
						{
							if ( extendedOverlaps[k].readStart > 0 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							break ;
						}
					}
					if ( j < i )
						continue ;

					// The previous overlaps might just failed to extend at one side.
					for ( j = 0 ; j < failedExtendedOverlapsCnt ; ++j )
					{
						if ( extendedOverlaps[k].readStart >= failedExtendedOverlaps[j].readStart &&
							extendedOverlaps[k].readEnd <= failedExtendedOverlaps[j].readEnd )
						{
							if ( extendedOverlaps[k].readStart > 0 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							if ( extendedOverlaps[k].readEnd < len - 1 )
							{
								if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor < 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
									seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor = 
										extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 ;
							}
							break ;
						}
					}
					if ( j < failedExtendedOverlapsCnt  )
						continue ;
					
					if ( extendedOverlaps[k].readStart > 0 )
					{
						if ( seqs[ extendedOverlaps[k].seqIdx ].minLeftExtAnchor >= extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
							continue ;
					}
					if ( extendedOverlaps[k].readEnd < len - 1 )
					{
						if ( seqs[ extendedOverlaps[k].seqIdx ].minRightExtAnchor >= extendedOverlaps[k].readEnd - extendedOverlaps[k].readStart + 1 )
							continue ;
					}


					tag = i ;
					++k ;
				}
				else
				{
					failedExtendedOverlaps[ failedExtendedOverlapsCnt ] = extendedOverlaps[k] ;
					++failedExtendedOverlapsCnt ;
				}
			}
			
			if ( k == 1 && 
				extendedOverlaps[0].readStart <= radius && extendedOverlaps[0].readEnd >= len - radius )
			{
				// Further check whehter this could merge two assemblies.
				// This happens when two contigs already overlap with each other.
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( tag == i )
						continue ;

					struct _seqWrapper &seq = seqs[ overlaps[i].seqIdx ] ; 
					if ( seq.isRef )
						continue ;
				
					if ( ExtendOverlap( r, len, seq, align, overlaps[i], extendedOverlaps[k] ) == 1 )
						++k ;
				}

				if ( k > 2 )
				{
					k = 1 ;
				}
				else if ( k == 2 )
				{
					if ( extendedOverlaps[0].seqEnd == seqs[ extendedOverlaps[0].seqIdx ].consensusLen - 1  
						&& extendedOverlaps[1].seqStart == 0  )
					{
						// no need to change.
						sortExtendedOverlaps = false ;
					}
					else if ( extendedOverlaps[0].seqStart == 0 && 
						extendedOverlaps[1].seqEnd == seqs[ extendedOverlaps[1].seqIdx ].consensusLen - 1 )
					{
						// swap 0, 1
						sortExtendedOverlaps = false ;
						
						struct _overlap tmp = extendedOverlaps[0] ; 
						extendedOverlaps[0] = extendedOverlaps[1] ;
						extendedOverlaps[1] = tmp ;
					}
					else
						k = 1 ;
				}
			}

			if ( similarityThreshold > novelSeqSimilarity )
			{
				// Filter extended sequence.
				int cnt = 0 ;
				for ( i = 0 ; i < k ; ++i )
				{
					if ( extendedOverlaps[i].similarity >= similarityThreshold )
					{
						extendedOverlaps[ cnt ] = extendedOverlaps[i] ;
						++cnt ;
					}
				}
				k = cnt ;
			}

#ifdef DEBUG
			for ( i = 0 ; i < k ; ++i )
				printf( "extended %d: %d %s. %d. %d %d %d %d %lf\n", i, extendedOverlaps[i].seqIdx, 
						seqs[ extendedOverlaps[i].seqIdx ].name, extendedOverlaps[i].strand, 
						extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, extendedOverlaps[i].seqStart, 
						extendedOverlaps[i].seqEnd, extendedOverlaps[i].similarity ) ; 
			fflush( stdout ) ;
#endif
			

			if ( k > 1 )
			{
				// If we are merging multiple novel seqs, make sure they are different seqs.
				for ( i = 0 ; i < k - 1 ; ++i )
					for ( j = i + 1 ; j < k ; ++j )
					{
						if ( extendedOverlaps[i].seqIdx == extendedOverlaps[j].seqIdx )
						{
							k = 0 ;
							break ;
						}
					}
			}

			/*if ( k == 1 ) 
			{
				// Check wether a match to the annotation fits better, in that case, creating 
				// 	a new sequence is a better.
				for ( i = 0 ; i < overlapCnt ; ++i )
				{
					if ( !seqs[ overlaps[i].seqIdx ].isRef )
						continue ;
					if ( overlaps[i].readStart == 0 && overlaps[i].readEnd == len - 1 
						&& overlaps[i].readStart <= extendedOverlaps[0].readStart 
						&& overlaps[i].readEnd >= extendedOverlaps[0].readEnd 
						&& overlaps[i].similarity > extendedOverlaps[0].similarity 
						&& overlaps[i].matchCnt > extendedOverlaps[0].matchCnt )
					{
						k = 0 ;
					}
				}
			}*/

			if ( k > 1 )
			{
				int eOverlapCnt = k ;
				addNew = false ;		
				// Merge sequences.
				// Reorder the overlaps to the order on the read coordinate.

				if ( sortExtendedOverlaps )
					std::sort( extendedOverlaps, extendedOverlaps + eOverlapCnt, CompSortOverlapsOnReadCoord ) ;
				
#ifdef DEBUG
				for ( i = 0 ; i < k ; ++i )
					printf( "sort extended %d: %d %s. %d. %d %d %d %d\n", i, extendedOverlaps[i].seqIdx, 
							seqs[ extendedOverlaps[i].seqIdx ].name, extendedOverlaps[i].strand, 
							extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, extendedOverlaps[i].seqStart, 
							extendedOverlaps[i].seqEnd ) ; 
#endif			
				// Infer from the gene name to make sure we would not merge two different receptors 
				for ( i = 0 ; i < k ; ++i )
				{
					for ( j = i + 1 ; j < k ; ++j )
					{
						if ( !IsNameCompatible( seqs[ extendedOverlaps[i].seqIdx ].name,
							seqs[ extendedOverlaps[j].seqIdx ].name ) )
						{
							delete[] extendedOverlaps ;
							delete[] failedExtendedOverlaps ;
							return -1 ;
						}
					}
				}
				// Compute the new consensus.
				int sum = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					sum += seqs[ extendedOverlaps[i].seqIdx ].consensusLen ;
				char *newConsensus = (char *)malloc( sizeof(char) * ( sum + len + 1 ) ) ;
				
				// Compute the location of the seqs in the new merged seq
				int *seqOffset = new int[ eOverlapCnt ] ; 
				int base = 0 ;
				
				if ( extendedOverlaps[0].readStart > 0 )
				{
					for ( i = 0 ; i < eOverlapCnt ; ++i )
						seqOffset[i] = extendedOverlaps[i].readStart ;
				}
				else
				{
					seqOffset[0] = 0 ;
					
					for ( i = 1 ; i < eOverlapCnt ; ++i )
					{
						seqOffset[i] = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - 1  
							+ ( extendedOverlaps[i].readStart - extendedOverlaps[i - 1].readEnd ) ;	
					}
				}
				
#ifdef DEBUG
				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					printf( "merge %d: %d %d %d %d %d. %d\n", i, extendedOverlaps[i].readStart, extendedOverlaps[i].readEnd, 
						extendedOverlaps[i].seqStart, extendedOverlaps[i].seqEnd, seqs[ extendedOverlaps[i].seqIdx ].consensusLen, 
						seqOffset[i] ) ;
				}
#endif

				// Copy the original consensus in.
				// The earlier seq has higher weight.
				for ( i = eOverlapCnt - 1 ; i >= 0 ; --i )
				{
					memcpy( newConsensus + seqOffset[i], seqs[ extendedOverlaps[i].seqIdx ].consensus,
						seqs[ extendedOverlaps[i].seqIdx ].consensusLen ) ;	
				}

				// Fill in the gaps
				if ( extendedOverlaps[0].readStart > 0 )
					memcpy( newConsensus, r, extendedOverlaps[0].readStart ) ;		
				
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( extendedOverlaps[i].readStart > extendedOverlaps[i - 1].readEnd + 1 )
						memcpy( newConsensus + seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen,
							r + extendedOverlaps[i - 1].readEnd + 1, 
							extendedOverlaps[i].readStart - extendedOverlaps[i - 1].readEnd - 1 ) ;
				}
				

				int newConsensusLen = 0 ;
				if ( extendedOverlaps[i - 1].readEnd < len )
				{
					memcpy( newConsensus + seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen,
							r + extendedOverlaps[i - 1].readEnd + 1,
							len - extendedOverlaps[i - 1].readEnd -1 ) ;
					newConsensusLen = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
					                     + ( len - extendedOverlaps[i - 1].readEnd - 1 ) ;
					newConsensus[ newConsensusLen ] = '\0' ;
				}
				else
				{
					newConsensusLen = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen ;
					newConsensus[ newConsensusLen ] = '\0' ;
				}
				//printf( "newconsensus %s %d\n", newConsensus, newConsensusLen ) ;
				
				// Rearrange the memory structure for posWeight.	
				int newSeqIdx = extendedOverlaps[0].seqIdx ;
				k = 0 ;
				for ( i = 1 ; i < eOverlapCnt ; ++i )
					if ( extendedOverlaps[i].seqIdx < newSeqIdx )
					{
						newSeqIdx = extendedOverlaps[i].seqIdx ;
						k = i ;
					}
				SimpleVector<struct _posWeight> &posWeight = seqs[newSeqIdx].posWeight ;
				posWeight.ShiftRight( seqOffset[k] ) ;
				posWeight.ExpandTo( newConsensusLen ) ;
				posWeight.SetZero( 0, seqOffset[k] ) ;
				posWeight.SetZero( seqOffset[k] + seqs[newSeqIdx].consensusLen, 
					newConsensusLen - seqs[ newSeqIdx ].consensusLen -seqOffset[k] ) ;
				

				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					// Though not the most efficient implementation, it seems very straightforward.
					int seqIdx = extendedOverlaps[i].seqIdx ;
					if ( seqIdx == newSeqIdx )
						continue ;

					for ( j = 0 ; j < seqs[ seqIdx ].consensusLen ; ++j )
					{
						int l ;
						posWeight[ seqOffset[i] + j ] += seqs[ seqIdx ].posWeight[j] ;
					}
				}

				// Update the index.
				KmerCode kmerCode( kmerLength ) ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					seqIndex.RemoveIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus,
						seqs[ extendedOverlaps[i].seqIdx ].consensusLen, extendedOverlaps[i].seqIdx, 0 ) ;
				
				/*if ( seqOffset[0] != 0 )
				{
					seqIndex.UpdateIndexFromRead( kmerCode, seqs[ extendedOverlaps[0].seqIdx ].consensus, 
							seqs[ extendedOverlaps[0].seqIdx].consensusLen, seqOffset[0], 
							extendedOverlaps[0].seqIdx, newSeqIdx ) ; 
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					//printf( "%d->%d.\n", extendedOverlaps[i].seqIdx, newSeqIdx ) ;
					// Don't use the overlapped portion, which will create duplicated entries in the index.
					int start = 0 ;
					if ( seqOffset[i] < seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - kmerLength + 1 )
					{
						start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen - kmerLength + 1 -
								seqOffset[i] ;
						printf( "start=%d\n", start ) ;
						seqIndex.RemoveIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus, 
								start + kmerLength - 1, extendedOverlaps[i].seqIdx, 0 ) ;
					}
					seqIndex.UpdateIndexFromRead( kmerCode, seqs[ extendedOverlaps[i].seqIdx ].consensus, 
							seqs[ extendedOverlaps[i].seqIdx].consensusLen, seqOffset[i], 
							extendedOverlaps[i].seqIdx, newSeqIdx ) ;
				}
				
				// Update the index for the gaps.
				if ( extendedOverlaps[0].readStart > 0 )
				{
					seqIndex.BuildIndexFromRead( kmerCode, r, extendedOverlaps[0].readStart + kmerLength - 1, newSeqIdx ) ;
				}
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( extendedOverlaps[i].readStart > extendedOverlaps[i - 1].readEnd - kmerLength + 1 )
					{
						int start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
									- kmerLength + 1 ;
						int end = seqOffset[i] + kmerLength - 2 ;
						int rstart = extendedOverlaps[i - 1].readEnd - kmerLength + 2 ;
						seqIndex.BuildIndexFromRead( kmerCode, r + rstart, end - start + 1, newSeqIdx, start ) ;
					}
				}
				if ( extendedOverlaps[i - 1].readEnd < len - 1 )
				{
					int start = seqOffset[i - 1] + seqs[ extendedOverlaps[i - 1].seqIdx ].consensusLen 
						- kmerLength + 1 ;
					int rstart = extendedOverlaps[i - 1].readEnd - kmerLength + 2 ;
					seqIndex.BuildIndexFromRead( kmerCode, r + rstart, len - rstart, newSeqIdx, start ) ;
				}
				*/

				// Update the name.
				// TODO: use array of names.
				sum = 0 ;
				for ( i = 0 ; i < eOverlapCnt ; ++i )
					sum += strlen( seqs[ extendedOverlaps[i].seqIdx ].name ) ;
				char* nameBuffer = new char[sum + eOverlapCnt + 1 ] ;
				
				strcpy( nameBuffer, seqs[ extendedOverlaps[0].seqIdx ].name ) ;
				sum = strlen( nameBuffer ) ;
				//printf( "%d\n", seqs.size() ) ;
				for ( i = 1 ; i < eOverlapCnt ; ++i )
				{
					if ( strcmp( seqs[ extendedOverlaps[i].seqIdx ].name, seqs[ extendedOverlaps[i - 1].seqIdx ].name ) )
					{
						nameBuffer[ sum ] = '+' ;
						nameBuffer[ sum + 1 ] = '\0' ;
						strcpy( nameBuffer + sum + 1, seqs[ extendedOverlaps[i].seqIdx ].name ) ;
						sum = sum + 1 + strlen( seqs[ extendedOverlaps[i].seqIdx ].name ) ;
					}
				}
				free( seqs[ newSeqIdx ].name ) ;
				seqs[ newSeqIdx ].name = strdup( nameBuffer ) ;
				delete[] nameBuffer ;

				// Relase the memory for merged seqs.
				for ( i = 0 ; i < eOverlapCnt ; ++i )
				{
					int seqIdx = extendedOverlaps[i].seqIdx ;
					if ( seqIdx == newSeqIdx )
						continue ;
					
					ReleaseSeq( seqIdx ) ;
				}
					
				// Put the new allocated stuff in.
				free( seqs[ newSeqIdx ].consensus ) ;
				seqs[ newSeqIdx ].consensus = newConsensus ;
				seqs[ newSeqIdx ].consensusLen = newConsensusLen ;
				
				// Update the index
				UpdateConsensus( newSeqIdx, false ) ;
				seqIndex.BuildIndexFromRead( kmerCode, newConsensus, newConsensusLen, newSeqIdx ) ;
				
				// Update the anchor requirement.
				seqs[ newSeqIdx ].minLeftExtAnchor = seqs[ extendedOverlaps[0].seqIdx ].minLeftExtAnchor ;
				seqs[ newSeqIdx ].minRightExtAnchor = seqs[ extendedOverlaps[ eOverlapCnt - 1 ].seqIdx ].minRightExtAnchor ;
				
				// either one of the ends of read or seq should be 0.
				readInConsensusOffset = 0 ;
				if ( extendedOverlaps[0].seqStart > 0 )
					readInConsensusOffset = extendedOverlaps[0].seqStart ;

				delete[] seqOffset ;

				seqIdx = newSeqIdx ;
			}
			else if ( k == 1 )
			{
				// Extend a sequence
				addNew = false ;

				seqIdx = extendedOverlaps[0].seqIdx ;
				struct _seqWrapper &seq = seqs[ extendedOverlaps[0].seqIdx ] ;
				
				// Compute the new consensus.
				if ( extendedOverlaps[0].readStart > 0 || extendedOverlaps[0].readEnd < len - 1 )
				{
					char *newConsensus = (char *)malloc( sizeof( char ) * ( 
						( extendedOverlaps[0].readStart + len - 1 -extendedOverlaps[0].readEnd ) + seq.consensusLen + 1 ) ) ;

					if ( extendedOverlaps[0].readStart > 0 )
					{
						// add read[0...readStart-1] to the consensus.
						//for ( i = 0 ; i < extendedOverlaps[0].readStart ; ++i )
						//	newConsensus[i] = r[i] ;
						memcpy( newConsensus, r, extendedOverlaps[0].readStart ) ;
					}
					memcpy( newConsensus + extendedOverlaps[0].readStart, seq.consensus, seq.consensusLen ) ;
					j = extendedOverlaps[0].readStart + seq.consensusLen ;	

					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						for ( i = extendedOverlaps[0].readEnd + 1 ; i < len ; ++i, ++j )
							newConsensus[j] = r[i] ;
					}
					newConsensus[j] = '\0' ;
					int newConsensusLen = strlen( newConsensus ) ;	
					//printf( "new consensus %s %d. %s %d\n", seq.consensus, seq.consensusLen, newConsensus, j ) ;
					
					// Update index 
					int shift = extendedOverlaps[0].readStart ;
					KmerCode kmerCode( kmerLength ) ;
					if ( shift > 0 )
					{
						seqIndex.BuildIndexFromRead( kmerCode, newConsensus, extendedOverlaps[0].readStart + kmerLength - 1, seqIdx ) ;
						seqIndex.UpdateIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, shift, seqIdx, seqIdx ) ; 
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int start = extendedOverlaps[0].readStart + extendedOverlaps[0].seqEnd - kmerLength + 2 ;
						seqIndex.BuildIndexFromRead( kmerCode, newConsensus + start , 
							( newConsensusLen - start ), seqIdx, 
							start ) ;
					}
					
					// Rearrange the memory structure for posWeight.
					int expandSize = extendedOverlaps[0].readStart + ( len - 1 - extendedOverlaps[0].readEnd ) ;
					seq.posWeight.ExpandBy( expandSize ) ;
					if ( shift > 0 )
					{
						for ( i = seq.consensusLen - 1 ; i >= 0 ; --i )
							seq.posWeight[i + shift] = seq.posWeight[i] ;
						
						seq.posWeight.SetZero( 0, shift ) ;
					}
					if ( extendedOverlaps[0].readEnd < len - 1 )
					{
						int start = extendedOverlaps[0].readStart + seq.consensusLen ;
						seq.posWeight.SetZero( start, len - extendedOverlaps[0].readEnd - 1 ) ;
					}
					
					// Update the anchor requirement
					if ( shift > 0 )
						seq.minLeftExtAnchor = 0 ;
					if ( extendedOverlaps[0].readEnd < len - 1 )
						seq.minRightExtAnchor = 0 ;

					// Adjust the name.
					// Find the possible ref seq
					int refIdx = -1 ;
					for ( i = 0 ; i < overlapCnt ; ++i )
					{
						if ( !seqs[ overlaps[i].seqIdx ].isRef )
							continue ;
						// Use refIdx as the idx in the overlaps list.
						if ( refIdx == -1 || 
							( overlaps[i].readEnd - overlaps[i].readStart > 
								overlaps[refIdx].readEnd - overlaps[refIdx].readStart ) )
						{
							refIdx = i ;
						}

						if ( strstr( seq.name, seqs[ overlaps[i].seqIdx ].name ) != NULL )
						{
							refIdx = i ;
							break ;
						}
					}
					if ( refIdx != -1 )
					{
						// Use refIdx as the idx in the seqs 
						refIdx = overlaps[ refIdx ].seqIdx ;
						if ( strstr( seq.name, seqs[ refIdx ].name ) == NULL )
						{
							char *nameBuffer = new char[ strlen( seqs[ refIdx ].name) + strlen( seq.name ) + 2 ] ;
							if ( extendedOverlaps[0].readStart > 0 )
							{
								sprintf( nameBuffer, "%s+%s", seqs[refIdx].name, seq.name ) ;
							}
							else
							{
								sprintf( nameBuffer, "%s+%s", seq.name, seqs[refIdx].name ) ;
							}
							free( seq.name ) ;
							seq.name = strdup( nameBuffer ) ;
							delete[] nameBuffer ;
						}
					}
					// either one of the ends of read or seq should be 0.
					readInConsensusOffset = 0 ;
					if ( extendedOverlaps[0].seqStart > 0 )
						readInConsensusOffset = extendedOverlaps[0].seqStart ;

					free( seq.consensus ) ;
					seq.consensus = newConsensus ;
					seq.consensusLen = newConsensusLen ;	
					//printf( "new consensus len %d\n", seq.consensusLen ) ;
				}
				else // the read is inside of the seq.
					readInConsensusOffset = extendedOverlaps[0].seqStart ;
			}

			// Update the posweight, assume we already compute the new readStart and shift existing posWeight.
			// seqIdx holds the index that we need to update.
			if ( !addNew )
			{
				struct _seqWrapper &seq = seqs[seqIdx] ;
				//printf( "%d %d. %d %d\n%s\n%s\n", seqIdx, seq.posWeight.Size(), readInConsensusOffset, len, seq.consensus, r) ;
				SimpleVector<int> nPos ;
				for ( i = 0 ; i < len ; ++i )
				{
					if ( r[i] == 'N' )
						continue ;
					++seq.posWeight[i + readInConsensusOffset].count[ nucToNum[ r[i] - 'A' ] ] ;
			
					if ( seq.consensus[i + readInConsensusOffset ] == 'N' )
					{
						nPos.PushBack( i ) ;
					}
				}
				SetPrevAddInfo( seqIdx, 0, len - 1, readInConsensusOffset, readInConsensusOffset + len - 1, overlaps[0].strand ) ;

				int size = nPos.Size() ;
				for ( i = 0 ; i < size ;  )
				{
					for ( j = i + 1 ; j < size ; ++j )
						if ( nPos[j] > nPos[j - 1] + kmerLength - 1 )
							break ;

					// [i,j) holds the N positions that are with kmer-length size.
					int l ;
					// Update the consensus
					for ( l = i ; l < j ; ++l )
						seq.consensus[ nPos[l] + readInConsensusOffset ] = r[ nPos[l] ] ;
					// Update the index
					KmerCode kmerCode( kmerLength ) ;
					int start = nPos[i] - kmerLength + 1 + readInConsensusOffset ;
					if ( start < 0 )
						start = 0 ;
					int end = nPos[j - 1] + kmerLength - 1 + readInConsensusOffset ;
					if ( end >= seq.consensusLen )
						end = seq.consensusLen - 1 ;
					seqIndex.BuildIndexFromRead( kmerCode, seq.consensus + start, end - start + 1, seqIdx, start  ) ;
					i = j ;
				}

				ret = seqIdx ;
			}
			
			free( rcRead ) ;
			delete[] align ;
		}

		k = 0 ;
		// See whether there is a reference seq match is sequence if it does not match any novel seq.
		int refSeqIdx = -1 ;
		if ( addNew )
		{
			for ( i = 0 ; i < overlapCnt ; ++i )
			{
				if ( seqs[ overlaps[i].seqIdx ].isRef )
				{
					refSeqIdx = overlaps[i].seqIdx ;
					break ;
				}
			}
			if ( i >= overlapCnt )
				addNew = false ;
			
			/*if ( !addNew )
			{
				for ( i = 0 ; i < overlapCnt ; ++i )
					if ( overlaps[i].similarity == 1.0 )
					{
						addNew = true ;
						refSeqIdx = -1 ;
						break ;
					}
			}*/
		}

		if ( addNew )
		{
			// A novel sequence
			// Go through the reference to annotate this read.
			for ( i = 0 ; i < overlapCnt ; ++i )					
			{
				// Check whether this overlap is used.
				for ( j = 0 ; j < k ; ++j )
				{
					;
				}
			}
			// TODO: If overlaps on C-gene, it should not be truncated.  
			
			// Add the sequence to SeqSet
			int idx = seqs.size() ;
			struct _seqWrapper ns ;
			
			if ( refSeqIdx >=0 )
				ns.name = strdup( seqs[ refSeqIdx ].name ) ;
			else
				ns.name = strdup( "unknown" ) ;
			ns.consensus = strdup( read ) ;
			ns.consensusLen = strlen( read ) ;
			if ( overlaps[0].strand == -1 )
				ReverseComplement( ns.consensus, read, len ) ;
			ns.isRef = false ;
			
			ns.posWeight.Reserve( len ) ;
			ns.posWeight.ExpandTo( len ) ;
			ns.posWeight.SetZero( 0, len ) ;
			for ( i = 0 ; i < len ; ++i )
			{
				//memset( ns.posWeight[i].count, 0, sizeof( ns.posWeight[i].count ) ) ;
				if ( ns.consensus[i] == 'N' )
					continue ;
				
				++ns.posWeight[i].count[ nucToNum[ ns.consensus[i] - 'A' ] ] ;
			}
			//printf( "%d %s %lld\n", ns.posWeight.Size(), ns.consensus, ns.posWeight.BeginAddress() ) ;
			ns.minLeftExtAnchor = ns.minRightExtAnchor = 0 ;
			seqs.push_back( ns ) ;

			// Don't forget to update index.
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, len, idx ) ;			
			
			SetPrevAddInfo( idx, 0, len - 1, 0, len - 1, overlaps[0].strand ) ; 
#ifdef DEBUG
			printf( "add novel seq: %d\n", idx ) ;
#endif
			ret = idx ;
		}

		delete[] extendedOverlaps ;
		delete[] failedExtendedOverlaps ;

		if ( ret == -1 )
		{
			SetPrevAddInfo( -2, -1, -1, -1, -1, 0 ) ; 
			ret = -2 ;
		}
		return ret ;
	}
	
	// Called when we just want to duplicate the add operation 
	//   applied before.
	int RepeatAddRead( char *read )
	{
		if ( prevAddInfo.seqIdx < 0 )
			return prevAddInfo.seqIdx ;

		int i ;
		char *r ;
				
		r = read ;
		if ( prevAddInfo.strand == -1 )
		{
			int len = strlen( read ) ;
			r = strdup( read ) ;
			ReverseComplement( r, read, len ) ;
		}
		
		struct _seqWrapper &seq = seqs[ prevAddInfo.seqIdx ] ;
		//printf( "%d %d. %d %d\n%s\n%s\n", seqIdx, seq.posWeight.Size(), readInConsensusOffset, len, seq.consensus, r) ;
		for ( i = prevAddInfo.readStart ; i <= prevAddInfo.readEnd ; ++i )
		{
			if ( r[i] == 'N' )
				continue ;
			++seq.posWeight[i + prevAddInfo.seqStart].count[ nucToNum[ r[i] - 'A' ] ] ;
		}

		if ( prevAddInfo.strand == -1 )
			free( r ) ;

		return prevAddInfo.seqIdx ;
	}

	void AddAssignedRead( char *read, struct _overlap assign )
	{
		if ( assign.seqIdx == -1 )
			return ;

		char *r = strdup( read ) ;
		int len = strlen( read ) ;
		if ( assign.strand == -1 )
			ReverseComplement( r, read, len ) ;
		
		UpdatePosWeightFromRead( seqs[ assign.seqIdx ].posWeight, assign.seqStart, r ) ;
		free( r ) ;
	}
	

	void UpdateAllConsensus()
	{
		int i ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
			UpdateConsensus( i, true ) ;
		}
	}

	void UpdateConsensus( int seqIdx, bool updateIndex )
	{
		int i, j ;
		struct _seqWrapper &seq = seqs[ seqIdx ] ;
		SimpleVector<struct _pair> changes ;
		for ( i = 0 ; i < seq.consensusLen ; ++i )
		{
			int max = 0 ;
			int maxTag = 0 ;
			for ( j = 0 ; j < 4 ; ++j )
				if ( seq.posWeight[i].count[j] > max )
				{
					max = seq.posWeight[i].count[j] ;
					maxTag = j ;
				}

			if ( max == 0 ) // A case of N
				continue ;

			if ( nucToNum[ seq.consensus[i] - 'A' ] != maxTag )
			{
				struct _pair np ;
				np.a = i ;
				np.b = maxTag ;
				changes.PushBack( np ) ;
			}
		}

		if ( changes.Size() == 0 )
			return ;
		
		if ( updateIndex )
		{
			// Inefficient implementation, improve in future.
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.RemoveIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, seqIdx, 0 ) ;
		}

		int size = changes.Size() ;
		for ( i = 0 ; i < size ; ++i )
			seq.consensus[ changes[i].a ] = numToNuc[ changes[i].b ] ;

		if ( updateIndex )
		{
			KmerCode kmerCode( kmerLength ) ;
			seqIndex.BuildIndexFromRead( kmerCode, seq.consensus, seq.consensusLen, seqIdx, 0 ) ;
		}
	}
	
	// Remove unneeded entries and rebuild the index.
	void Clean( bool removeRefSeq )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		k = 0 ;
		KmerCode kmerCode( kmerLength ) ;
		seqIndex.Clear() ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( seqs[i].consensus == NULL )
			{
				continue ;
			}
			if ( removeRefSeq && seqs[i].isRef )
			{
				//seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;
				
				ReleaseSeq( i ) ;
				continue ;
			}

			seqs[k] = seqs[i] ;
			//if ( i != k )
			//	seqIndex.UpdateIndexFromRead( kmerCode, seqs[k].consensus, seqs[k].consensusLen, 0, i, k ) ;
			seqIndex.BuildIndexFromRead( kmerCode, seqs[k].consensus, seqs[k].consensusLen, k, 0 ) ;
			++k ;
		}
		SetPrevAddInfo( -1, -1, -1, -1, -1, 0 ) ; 
		seqs.resize( k ) ;
	}

	void ChangeKmerLength( int kl )
	{
		kmerLength = kl ;
		Clean( false ) ;
	}
	
	// Find the seq id this read belongs to.
	int AssignRead( char *read, double similarity, struct _overlap &assign )
	{
		int i ;

		std::vector<struct _overlap> overlaps ;
		double backupSimilarity = novelSeqSimilarity ;
		novelSeqSimilarity = similarity ;
		
		int overlapCnt = GetOverlapsFromRead( read, 0, overlaps ) ;
		//printf( "%d %d\n", overlapCnt, mateOverlapCnt ) ;
		//printf( "%d %s\n%d %s\n", overlaps[0].strand, reads[i].seq, mateOverlaps[0].strand, reads[i + 1].seq ) ;
		assign.seqIdx = -1 ;

		if ( overlapCnt == 0 )
		{
			novelSeqSimilarity = backupSimilarity ;
			return -1 ;
		}
			
		std::sort( overlaps.begin(), overlaps.end() ) ;

		int len = strlen( read ) ;
		char *rc = new char[len + 1] ;

		ReverseComplement( rc, read, len ) ;

		char *r = read ;
		if ( overlaps[0].strand == -1 )
			r = rc ;

		struct _overlap extendedOverlap ;
		char *align = new char[ 2 * len + 2 ] ;
		/*for ( j = 0 ; j < overlapCnt ; ++j )
		  {
		  printf( "+ %d %d: %d %d %d %lf\n", i, j, overlaps[j].seqIdx, overlaps[j].seqStart, overlaps[j].seqEnd, overlaps[j].similarity) ;
		  }
		  for ( j = 0 ; j < mateOverlapCnt ; ++j )
		  {
		  printf( "- %d %d: %d %d %d %lf\n", i + 1, j, mateOverlaps[j].seqIdx, mateOverlaps[j].seqStart, mateOverlaps[j].seqEnd, mateOverlaps[j].similarity) ;
		  }*/
		int extendCnt = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			if ( ExtendOverlap( r, len, seqs[ overlaps[i].seqIdx ], align, 
						overlaps[i], extendedOverlap ) == 1 )
			{
				/*if ( extendCnt == 0 )
				  {
				  extendedOverlap = tmpExtendedOverlap ;
				  ++extendCnt ;
				  }
				  else if ( tmpExtendedOverlap.similarity == extend) */
				if ( extendedOverlap.readStart == 0 && extendedOverlap.readEnd == len - 1 )
					break ;
			}
		}

		novelSeqSimilarity = backupSimilarity ;
		if ( i >= overlapCnt )
		{
			delete[] rc ;
			delete[] align ;
			return -1 ;
		}
	
		delete[] rc ;
		delete[] align ;
		assign = extendedOverlap ;
		return assign.seqIdx ;
	}
	

	// Recompute the posweight based on assigned read
	void RecomputePosWeight( std::vector<struct _assignRead> &reads )
	{
		int i, j ;
		int readCnt = reads.size() ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
			seqs[i].posWeight.SetZero( 0, seqs[i].consensusLen ) ;

		for ( i = 0 ; i < readCnt ; ++i )
		{
			if ( reads[i].overlap.seqIdx == -1 )
				continue ;
			
			if ( reads[i].overlap.strand == 1 )
				UpdatePosWeightFromRead( seqs[ reads[i].overlap.seqIdx ].posWeight, reads[i].overlap.seqStart, reads[i].read ) ;
			else
			{
				char *r = strdup( reads[i].read ) ;
				ReverseComplement( r, reads[i].read, strlen( r ) ) ;
				UpdatePosWeightFromRead( seqs[ reads[i].overlap.seqIdx ].posWeight, reads[i].overlap.seqStart, r ) ;
				free( r ) ;
			}
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			struct _seqWrapper &seq = seqs[i] ;
			for ( j = 0 ; j < seq.consensusLen ; ++j )
				if ( seq.consensus[j] != 'N' && seq.posWeight[j].Sum() == 0 )
				{
					seq.posWeight[j].count[ nucToNum[ seq.consensus[j] - 'A' ] ] = 1 ;
				}
		}
	}
	
	// Return:the number of connections made.
	int ExtendSeqFromSeqOverlap( int leastOverlapLen )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		int ret = 0 ;
		repeatSimilarity = 0.99 ;

		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		struct _pair *next = new struct _pair[ seqCnt ] ; // a-index, b-the index in adj[i] 
		struct _pair *prev = new struct _pair[ seqCnt ] ;
		int *containedIn = new int[seqCnt] ;
	
		//BuildSeqOverlapGraph( 100, adj ) ;
		SimpleVector<bool> useInBranch ;
		useInBranch.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
			useInBranch[i] = true ;

		BuildBranchGraph( adj, leastOverlapLen, useInBranch ) ;
		// Keep only the connections that representing overlapping information.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].readEnd < seqs[i].consensusLen - 1 ||
					adj[i][j].seqStart > 0 )
					adj[i][j].seqIdx = -1 ;
				else if ( adj[i][j].readStart == 0 && adj[i][j].readEnd == seqs[i].consensusLen - 1
					&& adj[i][j].seqStart == 0 && adj[i][j].seqEnd == seqs[ adj[i][j].seqIdx ].consensusLen - 1 
					&& i > adj[i][j].seqIdx )
					adj[i][j].seqIdx = -1 ;
			}
			
			k = 0 ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].seqIdx != -1 )
				{
					adj[i][k] = adj[i][j] ;
					++k ;
				}
			}
			adj[i].resize( k ) ;
		}

		memset( next, -1, sizeof( struct _pair ) * seqCnt ) ;
		memset( prev, -1, sizeof( struct _pair ) * seqCnt ) ;
		memset( containedIn, -1, sizeof( int ) * seqCnt ) ;
		KmerCode kmerCode( kmerLength ) ;

		// Process the contained in.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			if ( size == 0 )
				continue ;
			
			if ( adj[i][0].readStart == 0 && adj[i][0].readEnd == seqs[i].consensusLen - 1 &&
				containedIn[ adj[i][0].seqIdx ] == -1 )
			{
				int seqIdx = adj[i][0].seqIdx ;
				//printf( "Contain: %d %d %lf\n%s\n%s\n", adj[i][0].readStart, adj[i][0].readEnd, adj[i][0].similarity, 
				//			seqs[i].consensus, seqs[ seqIdx ].consensus ) ;
				
				for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
				{
					seqs[ seqIdx ].posWeight[ adj[i][0].seqStart + j ] += seqs[i].posWeight[j] ;
				}

				seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;	
				ReleaseSeq( i ) ;

				containedIn[i] = seqIdx ;
				next[i].a = -2 ;
				prev[i].a = -2 ;
				
				++ret ;
			}
		}

		//Clean( true ) ;
		//return ret ;
		
		// Process the partial overlap case
		// Build the path
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			if ( size == 0 )
				continue ;
			k = 0 ;
			if ( containedIn[i] != -1 )
				continue ;

			if ( containedIn[ adj[i][0].seqIdx ] != -1 )
			{
				int father = containedIn[ adj[i][0].seqIdx ] ;
				while ( 1 )
				{
					if ( containedIn[ father ] == -1 )
						break ;
					father = containedIn[ father ] ;
				}
				for ( j = 0 ; j < size ; ++j )
					if ( adj[i][j].seqIdx == father )
						break ;
				if ( j < size )
					k = j ;
				else
					continue ;
			}
			if ( prev[ adj[i][k].seqIdx ].a == -1 )
			{
				next[i].a = adj[i][k].seqIdx ;
				next[i].b = k ;
				prev[ adj[i][k].seqIdx ].a = i ;
				prev[ adj[i][k].seqIdx ].b = k ;
			}
			else if ( prev[ adj[i][k].seqIdx ].a >= 0 )
			{
				next[ prev[ adj[i][k].seqIdx ].a ].a = -2 ;
				prev[ adj[i][k].seqIdx ].a = -2 ;
			}
		}
		
		//for ( i = 0 ; i < seqCnt ; ++i )
		//	printf( "%d next %d prev %d contained %d\n", i, next[i].a, prev[i].a, containedIn[i] ) ;
		
		// Use the paths to merge seqs.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( containedIn[i] == -1 && prev[i].a < 0 && next[i].a >= 0 ) // Head of a chain
			{
				std::vector<int> path ;
				int p = i ;
				int newConsensusLen = 0 ;
				std::vector<int> seqOffset ; 
				int offset = 0 ;

				k = 0 ;
				while ( 1 )
				{
					path.push_back( p ) ;
					++k ;
					seqOffset.push_back( offset ) ;
					
					if ( next[p].a >= 0 )
						offset += adj[p][ next[p].b ].readStart ;
					else
						break ;
					p = next[p].a ;
				}
			
				struct _seqWrapper ns ;
				ns.isRef = false ;

				// Obtain the length after merging.
				newConsensusLen = seqOffset[k - 1] + seqs[ path[k - 1] ].consensusLen ;
				char *newConsensus = (char *)malloc( sizeof( char ) * newConsensusLen + 1 ) ;
				for ( j = 0 ; j < k ; ++j )
					memcpy( newConsensus + seqOffset[j], seqs[ path[j] ].consensus, seqs[ path[j] ].consensusLen ) ;
				newConsensus[ newConsensusLen ] = '\0' ;
				
				ns.consensus = newConsensus ;
				ns.consensusLen = newConsensusLen ;
				ns.posWeight.ExpandTo( newConsensusLen ) ;
				
				// Update posweight
				ns.posWeight.SetZero( 0, newConsensusLen ) ;
				for ( j = 0 ; j < k ; ++j )
				{
					int l, c ;
					int seqIdx = path[j] ;
					for ( l = 0 ; l < seqs[ seqIdx ].consensusLen ; ++l )
						ns.posWeight[l + seqOffset[j] ] += seqs[ seqIdx ].posWeight[l] ; 
				} 
				// Update name
				int sum = 0 ;
				for ( j = 0 ; j < k ; ++j )
					sum += strlen( seqs[ path[j] ].name ) ;
				sum += k ;
				
				ns.name = (char *)malloc( sizeof( char ) * sum ) ;
				strcpy( ns.name, seqs[ path[0] ].name ) ;
				sum = strlen( ns.name ) ;
				for ( j = 1 ; j < k ; ++j )
				{
					ns.name[ sum ] = '+' ;
					ns.name[ sum + 1 ] = '\0' ;
					strcpy( ns.name + sum + 1, seqs[ path[j] ].name ) ;
					sum = sum + 1 + strlen( seqs[ path[j] ].name ) ;
				}

				// Update index	
				//int newSeqIdx = seqs.size() ;
				seqs.push_back( ns ) ;
				for ( j = 0 ; j < k ; ++j )
				{
					//seqIndex.RemoveIndexFromRead( kmerCode, seqs[ path[j] ].consensus, seqs[ path[j] ].consensusLen, path[j], 0 ) ;						
					ReleaseSeq( path[j] ) ;
				}
				//seqIndex.BuildIndexFromRead( kmerCode, ns.consensus, ns.consensusLen, newSeqIdx ) ;
			
				ret += k - 1 ;
			}
		}


		Clean( true ) ;
		
		delete[] adj ;
		delete[] next ;
		delete[] prev ;
		delete[] containedIn ;
		return ret ;
	}

	// Remove the seq that are substring of other seqs.
	int RemoveRedundantSeq()
	{
		int i ;
		int seqCnt = seqs.size() ;
		std::vector<struct _overlap> subsetOf ;
		subsetOf.resize( seqCnt ) ;
		
		for ( i = 0 ; i < seqCnt ; ++i )
			subsetOf[i].seqIdx = -1 ;	
		
		BuildSeqSubstringRelation( subsetOf ) ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( subsetOf[i].seqIdx != -1 )
				ReleaseSeq( i ) ;
		}
		Clean( true ) ;
	
		return seqs.size() ;
	}

	/*int AddRead( char *read )
	{
		int i, j, k ;
		int len = strlen( read ) ;

		std::vector<struct _overlap> overlaps ;
		int overlapCnt ;
		
		overlapCnt = GetOverlapsFromRead( read, overlaps ) ;
		
		if ( overlapCnt == 0 )
			return -1 ;
		
		struct _overlap *extendedOverlaps = new struct _overlap[ overlapCnt ] ;
		k = 0 ;
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			std::sort( overlaps.begin(), overlaps.end() ) ;
			if ( Extend)
		}
		int eOverlapCnt ;

		delete[] overlapCnt ;
	}*/

	void ResetPosWeight()
	{
		int i ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = seqs[i].posWeight.Size() ;
			seqs[i].posWeight.SetZero( 0, size ) ;
		}
	}
	
	int GetSeqCnt() 
	{
		return seqs.size() ;
	}

	int GetGeneType( char *name )
	{
		int geneType = -1 ;
		switch ( name[3] )
		{
			case 'V': geneType = 0 ; break ;
			case 'D': 
				  if ( name[4] >= '0' && name[4] <= '9' )
					  geneType = 1 ; 
				  else
					  geneType = 3 ;
				  break ;
			case 'J': geneType = 2 ; break ;
			default: geneType = 3 ; break ;
		}
		return geneType ;
	}

	bool IsSameGeneAllele( char *name, char *name2 )
	{
		int i ;
		int ret = true ;
		for ( i = 0 ; name[i] && name2[i] && name[i] != '*' && name2[i] != '*' ; ++i )
		{
			if ( name[i] != name2[i] )
			{
				ret = false ;
				break ;
			}
		}

		return ret ;
	}
	
	// The reference gene may have different length, which makes matchCnt criterion biased
	//   to longer gene, so we want to remove such effect
	// Return: is overlap a is better than b*threshold
	bool IsBetterGeneMatch( struct _overlap &a, struct _overlap &b, double threshold )
	{
		int matchCnt = a.matchCnt ;
		/*if ( a.seqEnd - a.seqStart + 1 > 0.95 * seqs[ a.seqIdx ].consensusLen ) 
			matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
		else if ( GetGeneType( seqs[ a.seqIdx ].name ) == 2 && a.seqEnd >= seqs[ a.seqIdx ].consensusLen - 3 
			&& a.readEnd >= seqs[a.seqIdx].consensusLen )
		{
			matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
		}*/
		int gapAllow = kmerLength + 1 ;
		if ( threshold >= 1 )
			gapAllow = 3 ;
		if ( a.seqIdx == -1 )
			return false ;
		if ( b.seqIdx == -1 )
			return true ;

		if ( GetGeneType( seqs[ a.seqIdx ].name ) == 2 /*&& threshold < 1*/ ) 
		{
			//printf( "hi %lf %lf %.1lf\n", a.similarity, b.similarity, threshold ) ;
			if ( a.seqEnd >= seqs[ a.seqIdx ].consensusLen - gapAllow /*&& a.readEnd >= seqs[a.seqIdx].consensusLen */
					&& b.seqEnd >= seqs[ b.seqIdx ].consensusLen - gapAllow /*&& b.readEnd >= seqs[b.seqIdx].consensusLen*/ )
			{
				//matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
				if ( a.similarity - 0.1 > b.similarity 
					&& a.seqEnd - a.seqStart > b.seqEnd - b.seqStart - kmerLength )
				{
					return true ;
				}
				else if ( a.similarity - 0.1 > b.similarity )
				{
					return true ;


					int adjustMatchCntA = matchCnt ;
					int adjustMatchCntB = b.matchCnt ;
					if ( a.readEnd - a.readStart < b.readEnd - b.readStart )
						adjustMatchCntA += ( b.readEnd - b.readStart - a.readEnd - a.readStart  ) / 2 ;
					else
						adjustMatchCntB += ( a.readEnd - a.readStart - b.readEnd - b.readStart  ) / 2 ;

					if ( adjustMatchCntA > adjustMatchCntB * threshold )
						return true ;
				}
				//else //if ( a.similarity > b.similarity )
				//	matchCnt = a.matchCnt / (double)seqs[ a.seqIdx ].consensusLen * seqs[ b.seqIdx ].consensusLen ;
			}
			else if (  a.seqEnd >= seqs[ a.seqIdx ].consensusLen - gapAllow && a.readEnd >= seqs[a.seqIdx].consensusLen ) 
			{	
				return true ;
			}
		}
		
		//if ( threshold == 1 )
		//	printf( "%s %d %lf %s %d %lf\n", seqs[ a.seqIdx ].name, a.matchCnt, a.similarity,
		//		 seqs[ b.seqIdx ].name, b.matchCnt, b.similarity ) ;

		if ( a.readStart == b.readStart && a.readEnd == b.readEnd )
		{
			if ( a.similarity > b.similarity )
				return true ;
			else if ( a.similarity < b.similarity )
				return false ;
		}

		if ( matchCnt > b.matchCnt * threshold )
			return true ;
		else if ( threshold < 1.0 && a.matchCnt + 10 >= b.matchCnt )
			return true ;
		else
			return false ;
	}
	
	// Figure out the gene composition for the read. 
	// Return successful or not.
	int AnnotateRead( char *read, int detailLevel, struct _overlap geneOverlap[4], struct _overlap cdr[3], char *buffer )
	{
		int i, j, k ;
		
		std::vector<struct _overlap> overlaps ;
		std::vector<struct _overlap> allOverlaps ;
		int overlapCnt ;
		
		char BT = '\0' ; // Bcell, Tcell
		char chain = '\0' ;
		int len = strlen( read ) ;
		buffer[0] = '\0' ;

		geneOverlap[0].seqIdx = geneOverlap[1].seqIdx = geneOverlap[2].seqIdx = geneOverlap[3].seqIdx = -1 ;
		if ( detailLevel >= 2 )
			cdr[0].seqIdx = cdr[1].seqIdx = cdr[2].seqIdx = -1 ;

		//sprintf( buffer, "%d", len ) ;
		if ( detailLevel > 0 )
			hitLenRequired = 17 ;	
		overlapCnt = GetOverlapsFromRead( read, detailLevel == 0 ? 0 : 1, overlaps ) ;		
		hitLenRequired = 31 ;

		if ( overlapCnt == 0 )
		{
			if ( detailLevel >= 2 )
				sprintf( buffer + strlen( buffer ), " * * * CDR1(0-0):0.00=null CDR2(0-0):0.00=null CDR3(0-0):0.00=null" ) ;
			return 0 ;
		}
		std::sort( overlaps.begin(), overlaps.end() ) ;
		// Get the coverage of the genes.
		for ( i = 0 ; i < overlapCnt ; ++i )
		{
			char *name = seqs[ overlaps[i].seqIdx ].name ;
			if ( BT && name[0] != BT )
				continue ;
			BT = name[0] ;
			
			if ( chain && name[2] != chain )
				continue ;
			chain = name[2] ;

			int geneType = GetGeneType( name ) ;
			//printf( "%d %s %lf %d: %d %d\n", i, seqs[ overlaps[i].seqIdx].name, overlaps[i].similarity, overlaps[i].matchCnt, overlaps[i].readStart, overlaps[i].readEnd ) ;
			
			if ( geneType >= 0 && geneOverlap[ geneType ].seqIdx == -1 )
				geneOverlap[ geneType ] = overlaps[i] ;
			
			//printf( "%s %d\n", seqs[ overlaps[i].seqIdx ].name, overlaps[i].matchCnt ) ;
			if ( geneType >= 0 && IsBetterGeneMatch( overlaps[i], geneOverlap[geneType], 0.95 ) )
			{
				allOverlaps.push_back( overlaps[i] ) ;
			}
			else if ( geneType >= 0 && geneOverlap[ geneType ].seqIdx != -1 &&
				( overlaps[i].readEnd < geneOverlap[geneType].readStart || 
					overlaps[i].readStart > geneOverlap[ geneType ].readEnd )
				&& IsBetterGeneMatch( overlaps[i], geneOverlap[geneType], 0.9 ) )
			{
				// The gene on a different region of the read. Might happen due to false alignment.
				allOverlaps.push_back( overlaps[i] ) ;
			}
		}
		
		// Check whether the match to the constant gene is random.
		if ( geneOverlap[3].seqIdx != -1 && geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 <= len / 2 
			&& geneOverlap[3].readEnd - geneOverlap[3].readStart + 1 <= 50 )
		{
			for ( i = 0 ; i < 3 ; ++i )
			{
				if ( geneOverlap[i].seqIdx >= 0 && 
					( geneOverlap[i].readEnd - 17 > geneOverlap[3].readStart 
						|| geneOverlap[3].readEnd < geneOverlap[i].readEnd ) 
					&& geneOverlap[3].seqStart >= 100 )
				{
					// Filter out all the overlaps from C gene
					geneOverlap[3].seqIdx = -1 ;
					break ;
				}
			}

			if ( i < 3 && detailLevel >= 1 )
			{
				int size = allOverlaps.size() ;
				for ( i = 0, k = 0 ; i < size; ++i )
				{
					int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
					if ( geneType != 3 )
					{
						allOverlaps[k] = allOverlaps[i] ;
						++k ;
					}
				}
				allOverlaps.resize( k ) ;
			}
		}
		
		// Extend overlap
		if ( detailLevel >= 1 )
		{
			char *align = new char[ 2 * len + 2 ] ;
			char *rvr = new char[len + 1] ;
			int size = allOverlaps.size() ;
			for ( i = 0 ; i < size ; ++i )
			{
				// Extend right.
				int seqIdx = allOverlaps[i].seqIdx ;				
				AlignAlgo::GlobalAlignment_OneEnd( seqs[ seqIdx ].consensus + allOverlaps[i].seqEnd + 1, 
					seqs[ seqIdx ].consensusLen - allOverlaps[i].seqEnd - 1, 
					read + allOverlaps[i].readEnd + 1, len - allOverlaps[i].readEnd - 1, 0, align ) ;
				//AlignAlgo::VisualizeAlignment( seqs[ seqIdx ].consensus + allOverlaps[i].seqEnd + 1, 
				//	seqs[ seqIdx ].consensusLen - allOverlaps[i].seqEnd - 1, 
				//	read + allOverlaps[i].readEnd + 1, len - allOverlaps[i].readEnd - 1, align ) ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						++allOverlaps[i].readEnd ;
						++allOverlaps[i].seqEnd ;
						
						if ( align[j] == EDIT_MATCH )
							allOverlaps[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							++allOverlaps[i].readEnd ;
						else if ( align[j] == EDIT_DELETE )
							++allOverlaps[i].seqEnd ;
					}
					else
						break ;
				}

				// Extend left.
				char *rvs = new char[seqs[ seqIdx ].consensusLen ] ;
				Reverse( rvr, read, allOverlaps[i].readStart ) ;
				Reverse( rvs, seqs[seqIdx].consensus, allOverlaps[i].seqStart ) ;
				//rvr[geneOverlap[i].readStart] = '\0' ;
				//rvs[geneOverlap[i].seqStart] = '\0' ;
				
				AlignAlgo::GlobalAlignment_OneEnd( rvs, allOverlaps[i].seqStart, rvr, allOverlaps[i].readStart, 0, align ) ;
				//AlignAlgo::VisualizeAlignment( rvs, geneOverlap[i].readStart, rvr, rvr[geneOverlap[i].readStart], align ) ;
				for ( j = 0 ; align[j] != -1 ; ++j )
				{
					if ( align[j] == EDIT_MATCH || align[j] == EDIT_MISMATCH )
					{
						--allOverlaps[i].readStart ;
						--allOverlaps[i].seqStart ;
					
						if ( align[j] == EDIT_MATCH )
							allOverlaps[i].matchCnt += 2 ;
					}
					else if ( radius > 0 )
					{
						if ( align[j] == EDIT_INSERT )
							--allOverlaps[i].readStart ;
						else if ( align[j] == EDIT_DELETE )
							--allOverlaps[i].seqStart ;
					}
					else
						break ;
				}
				delete[] rvs ;
				
				allOverlaps[i].similarity = (double)( allOverlaps[i].matchCnt ) / 
						( allOverlaps[i].seqEnd - allOverlaps[i].seqStart + 1 + 
							allOverlaps[i].readEnd - allOverlaps[i].readStart + 1 ) ;
			}
			
			for ( i = 0 ; i < 4 ; ++i )
			{
				geneOverlap[i].seqIdx = -1 ;
				geneOverlap[i].matchCnt = -1 ;
			}
			
			for ( i = 0 ; i < size ; ++i )
			{
				int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
				if ( IsBetterGeneMatch( allOverlaps[i], geneOverlap[ geneType ], 1.0 ) )
				{
					//if ( geneOverlap[ geneType].seqIdx != -1 )
					//	printf( "%s %s\n", seqs[ geneOverlap[ geneType ].seqIdx ].name, seqs[ allOverlaps[i].seqIdx ].name ) ;
					//else
					//	printf( "-1 %s\n", seqs[ allOverlaps[i].seqIdx ].name ) ;
						
					geneOverlap[ geneType ] = allOverlaps[i] ;
				}
			}

			// Rescue constant gene if the anchor is too short
			/*if ( geneOverlap[2].seqIdx != -1 && geneOverlap[3].seqIdx == -1 && 
				geneOverlap[2].readEnd + hitLenRequired - 1 >= len )
			{
				int seqCnt = seqs.size() ;
				int rstart = geneOverlap[2].readEnd + 6 ;
				int plen = len - rstart ;
				if ( plen >= 9 )
				{
					char *tBuffer = new char[len - geneOverlap[2].readEnd + 6] ;
					for ( j = 0 ; j < seqCnt ; ++j )
					{
						int geneType = GetGeneType( seqs[j].name ) ;
						if ( geneType != 3 )
							continue ;
						memcpy( tBuffer, seqs[j].consensus, len - geneOverlap[2].readEnd + 2 ) ;
						tBuffer[ len - geneOverlap[2].readEnd + 2 ] = '\0' ;
						char *p = strstr( tBuffer, read + rstart ) ;
						if ( p == NULL )
							continue ;
						int l ;
						for ( l = rstart ; l >= 0 && p >= tBuffer ; --p, --l )
							if ( read[l] != *p )
							{
								break ;
							}
						if ( p >= tBuffer )
							continue ;

						struct _overlap no ;
						no.seqIdx = j ;
						no.readStart = l ;
						no.readEnd = len - 1 ;
						no.seqStart = 0 ;
						no.seqEnd = no.readEnd - no.readStart ;
						no.matchCnt = 2 * ( no.readEnd - no.readStart ) ;
						no.similarity = 1.0 ;
						geneOverlap[ geneType ] = no ;
						allOverlaps.push_back( no ) ;
					}
					delete[] tBuffer ;
				}
			}*/
			
			// Test whether the V gene's coordinate is wrong if we have good J,C gene alignment.
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 && geneOverlap[3].seqIdx != -1 )
			{
				if ( geneOverlap[2].readEnd + 3 >= geneOverlap[3].readStart && 
					geneOverlap[2].readEnd - 3 <= geneOverlap[3].readStart &&
					geneOverlap[0].readEnd > geneOverlap[2].readStart + 6 )
				{
					struct _overlap orig = geneOverlap[0] ;
					geneOverlap[0].seqIdx = -1 ;
					geneOverlap[0].matchCnt = -1 ;
					for ( i = 0 ; i < size ; ++i )
					{
						int geneType = GetGeneType( seqs[ allOverlaps[i].seqIdx ].name ) ;
						if ( geneType != 0 )
							continue ;
						
						if ( allOverlaps[i].readEnd <= geneOverlap[2].readStart + 6 && 
							IsBetterGeneMatch( allOverlaps[i], geneOverlap[ geneType ], 1.0 ) )
						{
							geneOverlap[ geneType ] = allOverlaps[i] ;
						}
					}
					//if ( geneOverlap[0].seqIdx == -1 )
					//	geneOverlap[0] = orig ;
				}
			}

			delete[] align ;
			delete[] rvr ;
		}

		// Infer CDR1,2,3.
		char *cdr1 = NULL ;
		char *cdr2 = NULL ;
		
		if ( detailLevel >= 2 && geneOverlap[0].seqIdx != -1 
			&& ( geneOverlap[2].seqIdx == -1 || geneOverlap[0].readStart < geneOverlap[2].readStart ) )
		{
			// Infer CDR1, 2
			char *align = new char[ 2 * len + 2 ] ;
			struct _overlap vgene = geneOverlap[0] ; // Overlap with v-gene
			AlignAlgo::GlobalAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, vgene.seqEnd - vgene.seqStart + 1,
				read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, align ) ;
			//AlignAlgo::VisualizeAlignment( seqs[ vgene.seqIdx ].consensus + vgene.seqStart, vgene.seqEnd - vgene.seqStart + 1,
			//	read + vgene.readStart, vgene.readEnd - vgene.readStart + 1, align ) ;

			// Locate CDR1.
			int cdrIdx ;
			for ( cdrIdx = 0 ; cdrIdx <= 1 ; ++cdrIdx )
			{
				int seqRangeStart = seqs[ vgene.seqIdx ].info[ cdrIdx ].a ; 
				int seqRangeEnd = seqs[ vgene.seqIdx ].info[ cdrIdx ].b ;

				if ( vgene.seqStart <= seqRangeStart && vgene.seqEnd >= seqRangeEnd )
				{
					i = vgene.readStart - 1 ;
					j = vgene.seqStart - 1 ;
					int readRangeStart, readRangeEnd ;
					int matchCnt = 0 ;
					for ( k = 0 ; align[k] != -1 ; ++k )	
					{
						if ( align[k] != EDIT_DELETE )
							++i ;
						if ( align[k] != EDIT_INSERT )
							++j ;
						
						if ( j == seqRangeStart )
							readRangeStart = i ;
						if ( j >= seqRangeStart && align[k] == EDIT_MATCH )
							matchCnt += 2 ;
						if ( j == seqRangeEnd )
						{
							readRangeEnd = i ;
							break ;
						}
					}
					cdr[cdrIdx].seqIdx = vgene.seqIdx ;
					cdr[cdrIdx].readStart = readRangeStart ;
					cdr[cdrIdx].readEnd = readRangeEnd ;
					cdr[cdrIdx].matchCnt = matchCnt ;
					cdr[cdrIdx].similarity = (double)matchCnt / 
						( readRangeEnd - readRangeStart + 1 + seqRangeEnd - seqRangeStart + 1 ) ;
					//printf( "%d: %d %d; %d %d\n", matchCnt, readRangeStart, readRangeEnd, seqRangeStart, seqRangeEnd ) ;

					char *r =  ( char * )malloc( sizeof( char ) * ( readRangeEnd - readRangeStart + 2 ) ) ;
					memcpy( r, read + readRangeStart, readRangeEnd - readRangeStart + 1 ) ;
					r[  readRangeEnd - readRangeStart + 1 ] = '\0' ;
					if ( cdrIdx == 0 )
						cdr1 = r ;
					else if ( cdrIdx == 1 )
						cdr2 = r ;
				}
			}
			delete[] align ;
		}
		
		char *cdr3 = NULL ;
		double cdr3Score = 0 ;
		
		if ( detailLevel >= 2 )
		{
			// Infer CDR3.
			int s, e ;
			int boundS = 0, boundE = len - 2 ; //[boundS, boundE)
			int range = 37 ;
			if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 )
			{
				// The case that we have anchor.
				// Find the motif for anchor
				if ( geneOverlap[2].readEnd > geneOverlap[0].readEnd )
				{
					int startFrame = geneOverlap[0].seqStart % 3 ;
					int ns = geneOverlap[0].readEnd ; //+ ( seqs[ geneOverlap[0].seqIdx ].consensusLen - 1 - geneOverlap[0].seqEnd ) ;
					s = ns - ( ns - geneOverlap[0].readStart + startFrame ) % 3 ;
					s = s + 6 < len ? s + 6 : s ;
					startFrame = ( seqs[ geneOverlap[2].seqIdx ].consensusLen - 1 - geneOverlap[2].seqEnd ) % 3 ;
					//e = geneOverlap[2].readStart + 
					//	( geneOverlap[2].readEnd - geneOverlap[2].readStart + startFrame ) % 3 ;
					e = geneOverlap[2].readStart ;
					e = e - 6 >= 0 ? e - 6 : e ;
					int adjustE = e ;
					int locate = -1 ;
					for ( i = adjustE ; i < geneOverlap[2].readEnd  && i + 11 < len ; ++i )
					{
						// A strong motif, should be used as the anchor.
						if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
									DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
								&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locate = i ;
							break ;
						}
					}

					if ( locate != -1 )
						e = locate ;

					if ( e < s + 12 )
						range += 15 ;

					if ( s - range > boundS )
						boundS = s - range ;
					if ( e + range < boundE )
						boundE = e + range ;

					if ( locate != -1 )
						s = s + ( e - s ) % 3 ;
				}
				else
				{
					s = 0 ;
					e = len ;
					boundS = 1 ;
				}	
			
			}
			else if ( geneOverlap[2].seqIdx != -1 )
			{
				e = geneOverlap[2].readStart ;
				e = e - 6 >= 0 ? e - 6 : e ;
				int adjustE = e ;
				int locate = -1 ;
				for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
				{
					if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
								DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
							&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
					{
						locate = i ;
						break ;
					}
				}
				if ( locate != -1 )
					e = locate ;

				s = e - 12 ;
				if ( s < 0 )
					s = 0 ;
				if ( e + 31 < boundE )
					boundE = e + 31 ;
			}
			else if ( geneOverlap[0].seqIdx != -1 
				&& geneOverlap[0].seqEnd >= seqs[ geneOverlap[0].seqIdx ].consensusLen - 50 )
			{
				int startFrame = geneOverlap[0].seqStart % 3 ;
				s = geneOverlap[0].readEnd + ( geneOverlap[0].readEnd - geneOverlap[0].readStart - startFrame ) % 3 ;
				s = s + 6 < len ? s + 6 : s ;
				e = s + 12 ;
				if ( s - 31 > boundS )
					boundS = s - 31 ;
				
				int adjustE = e ;
				int locate = -1 ;
				if ( geneOverlap[3].seqIdx != -1 )
					boundE = geneOverlap[3].readStart - 2 ;
				for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
				{
					if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
								DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
							&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
					{
						locate = i ;
						break ;
					}
				}

				if ( locate != -1 )
				{
					e = locate ;
					s = e - 12 ;
					if ( s < 0 )
						s = 0 ;
				}
			}
			else
			{
				s = 0 ;
				e = len ;
				boundS = 1 ;
			}
			
			if ( geneOverlap[2].seqIdx != -1 && boundE > geneOverlap[2].readEnd )
				boundE = geneOverlap[2].readEnd ;
				

			int locateS = -1 ;
			int locateE = -1 ;
			
			// The YYC motif on V gene, mentioned in TRUST3 paper, but seems not mentioned in IMGT Junction Analysis.
			for ( i = s ; i >= boundS ; i -= 3 )
			{
				if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
					&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
					&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
				{
					locateS = i + 6 ;
					break ;
				}
			}
			

			if ( locateS == -1 )
			{
				// Don't follow the frame rule 
				for ( i = s ; i >= boundS ; --i )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
							&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
					{
						locateS = i + 6 ;
						break ;
					}
				}
			}
			
			if ( locateS == -1 )
			{
				// Try the YxC motif
				for ( i = s ; i >= boundS ; i -= 3 )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							//&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
							&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
					{
						locateS = i + 6 ;
						break ;
					}
				}
				
				if ( locateS == -1 )
				{
					for ( i = s ; i >= boundS ; --i )
					{
						if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
								//&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' 
								&& DnaToAa( read[i + 6], read[i + 7], read[i + 8] ) == 'C' )
						{
							locateS = i + 6 ;
							break ;
						}
					}

				}
			}

			if ( locateS == -1 )
			{
				for ( i = s ; i >= boundS ; i -= 3 )
				{
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
					{
						locateS = i ;
						break ;
					}
				}
			}
			
			
			if ( locateS == -1 && geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx != -1 )
			{
				// Expand the search range.
				int newS = e - 12 ; //- ( e - 12 - s ) % 3 ;
				if ( newS > s )
				{
					for ( i = newS ; i > s ; i -= 3 )
						if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
						{
							locateS = i ;
							break ;
						}
					if ( locateS == -1 )
					{
						for ( i = newS ; i > s ; --i )
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
							{
								locateS = i ;
								break ;
							}
					}

				}
			}

			if ( locateS == -1 )
			{
				for ( i = s ; i >= boundS ; --i )
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'C' ) 
					{
						locateS = i ;
						break ;
					}
			}

			if ( locateS == -1 )
			{
				// YYx motif.
				for ( i = s ; i >= boundS ; --i )
					if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'Y' 
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] )== 'Y' )
					{
						locateS = i + 6 ;
						break ;
					}
			}
			
			int adjustE = e  ;
			if ( 1 ) //locateS != -1 )
			{
				if ( locateS != -1 )
					adjustE = e - ( e - locateS ) % 3 ; 
				for ( i = adjustE ; i < boundE && i + 11 < len ; i += 3 )
				{
					if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
								DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
							&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
							&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
					{
						locateE = i ;
						break ;
					}
				}

				if ( locateE == -1 )
				{
					adjustE = e ;
					for ( i = adjustE ; i < boundE && i + 11 < len ; ++i )
					{
						if ( ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' || 
									DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )  
								&& DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&& DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locateE = i ;
							break ;
						}
					}
				}

				if ( locateE == -1 )
				{
					// Use weaker motif.
					if ( locateS != -1 )
					{
						adjustE = e - ( e - locateS ) % 3 ;
						if ( adjustE + 3 < locateS + 18 )
						adjustE = locateS + 15 ;
					}
					for ( i = adjustE ; i < boundE ; ++i )
					{
						// The GxG motif
						if ( read[i] == 'T' && DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
								&&  DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
						{
							locateE = i ;
							break ;
						}
					}


					if ( locateE == -1 )
					{
						for ( i = adjustE ; i < boundE ; ++i )
						{
							// The GxG motif
							if ( read[i] == 'T' && DnaToAa( read[i + 3], read[i + 4], read[i + 5] ) == 'G' 
									&&  DnaToAa( read[i + 9], read[i + 10], read[i + 11] ) == 'G' )
							{
								locateE = i ;
								break ;
							}
						}

					}
					
					if ( locateE == -1 )
					{
						for ( i = adjustE ; i < boundE ; i += 3 )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' )
							{
								//printf( "%c%c%c=>%c\n", read[i], read[i + 1], read[i + 2],
								//	DnaToAa( read[i], read[i + 1], read[i + 2] ) ) ;
								//printf( "%c%c%c=>%c\n", read[j], read[j + 1], read[j + 2],
								//	DnaToAa( read[j], read[j + 1], read[j + 2] ) ) ;
								locateE = i ;
								break ;
							}
						}
					}
					if ( locateE == -1 )
					{
						for ( i = e ; i < boundE ; i += 3 )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )
							{
								locateE = i ;
								break ;
							}
						}
					}

					if ( locateE == -1 )
					{
						// frame shift happens or no locateS.
						adjustE = e ; 
						for ( i = adjustE ; i < boundE ; ++i )
						{
							if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'W' )
							{
								//printf( "%c%c%c=>%c\n", read[i], read[i + 1], read[i + 2],
								//	DnaToAa( read[i], read[i + 1], read[i + 2] ) ) ;
								//printf( "%c%c%c=>%c\n", read[j], read[j + 1], read[j + 2],
								//	DnaToAa( read[j], read[j + 1], read[j + 2] ) ) ;
								locateE = i ;
								break ;
							}
						}
						if ( locateE == -1 )
						{
							for ( i = e ; i < boundE ; ++i )
							{
								if ( DnaToAa( read[i], read[i + 1], read[i + 2] ) == 'F' )
								{
									locateE = i ;
									break ;
								}
							}
						}
					}
				}
			}

			if ( locateE + 2 - locateS + 1 < 18 )
			{
				if ( geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 )
					locateS = -1 ;
				else if ( geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1  )
					locateE = -1 ;
			}

			if ( locateS != -1 && locateE != -1 && locateE + 2 - locateS + 1 >= 18 )
			{
				s = locateS ;
				e = locateE + 2 ;

				cdr3 = new char[e - s + 2  + 1 ] ;
				memcpy( cdr3, read + s, e - s + 1 ) ;
				cdr3[e - s + 1] = '\0' ;

				cdr[2].seqIdx = 0 ;
				cdr[2].readStart = s ;
				cdr[2].readEnd = e ;

				// Use the anchor motif to score the cdr
				if ( locateS - 6 > 0 )
					if ( DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) == 'Y' )
						cdr3Score += 100.0 / 6 ;
				if ( locateS - 3 > 0 )
					if ( DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) == 'Y' )
						cdr3Score += 100.0 / 6 ;

				if ( DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2 ] ) == 'C' )
					cdr3Score += 100.0 / 6 ;

				if ( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
					 DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' )
					cdr3Score += 100.0 / 6 ;

				if ( locateE + 5 < len )
					if ( DnaToAa( read[locateE + 3], read[ locateE + 4], read[ locateE + 5 ] ) == 'G' )
						cdr3Score += 100.0 / 6 ;
				if ( locateE + 11 < len )
					if ( DnaToAa( read[locateE + 9], read[ locateE + 10], read[ locateE + 11 ] ) == 'G' )
						cdr3Score += 100.0 / 6 ;

			}
			// Partial CDR3s
			else if ( locateS == -1 && locateE != -1 && geneOverlap[0].seqIdx == -1 && geneOverlap[2].seqIdx != -1 
				&& locateE + 11 < len && locateE > 15 && locateE <= 60 ) 
			{
				if ( ( DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'W' ||
					 DnaToAa( read[locateE], read[ locateE + 1], read[ locateE + 2 ] ) == 'F' )
					 && DnaToAa( read[locateE + 3], read[ locateE + 4], read[ locateE + 5 ] ) == 'G' 
					 && DnaToAa( read[locateE + 9], read[ locateE + 10], read[ locateE + 11 ] ) == 'G' )
				{
					locateS = locateE % 3 ;			
					cdr3Score = 0 ;

					s = locateS ;
					e = locateE + 2 ;
					
					if ( e - s + 1 >= 18 )
					{
						cdr3 = new char[e - s + 2  + 1 ] ;
						memcpy( cdr3, read + s, e - s + 1 ) ;
						cdr3[e - s + 1] = '\0' ;

						cdr[2].seqIdx = 0 ;
						cdr[2].readStart = s ;
						cdr[2].readEnd = e ;
					}
				}
			}
			else if ( locateS != -1 && locateE == -1 && geneOverlap[0].seqIdx != -1 && geneOverlap[2].seqIdx == -1 
				&& locateS - 6 > 0 && locateS + 18 < len && locateS + 2 + 60 > len )
			{
				if ( DnaToAa( read[locateS], read[ locateS + 1], read[ locateS + 2 ] ) == 'C' 
					&& DnaToAa( read[locateS - 3], read[ locateS - 2], read[ locateS - 1 ] ) == 'Y'  
					&& DnaToAa( read[locateS - 6], read[ locateS - 5], read[ locateS - 4 ] ) == 'Y' )
				{
					locateE = len - 3 - ( len - 3 - locateS ) % 3 ;
					cdr3Score = 0 ;
				
					s = locateS ;
					e = locateE + 2 ;
					
					if ( e - s + 1 >= 18 )
					{
						cdr3 = new char[e - s + 2  + 1 ] ;
						memcpy( cdr3, read + s, e - s + 1 ) ;
						cdr3[e - s + 1] = '\0' ;

						cdr[2].seqIdx = 0 ;
						cdr[2].readStart = s ;
						cdr[2].readEnd = e ;
					}
				}
			}
				
		}

		// Compute the name
		for ( i = 0 ; i < 4 ; ++i )
		{	
			if ( i == 1 ) // skip the D gene.
				continue ;
			if ( geneOverlap[i].seqIdx != -1 )
			{
				int offset = strlen( buffer ) ;
				int seqIdx = geneOverlap[i].seqIdx ;
				sprintf( buffer + offset, " %s(%d):(%d-%d):(%d-%d):%.2lf",
						seqs[ seqIdx ].name, seqs[ seqIdx ].consensusLen,
						geneOverlap[i].readStart, geneOverlap[i].readEnd, 
						geneOverlap[i].seqStart, geneOverlap[i].seqEnd, geneOverlap[i].similarity * 100 ) ;

				// Output the secondary assignment	
				int size = allOverlaps.size() ;
				int reportCnt = 0 ;
				SimpleVector<int> usedSeqIdx ;
				usedSeqIdx.Reserve( 5 ) ;
				for ( j = 0 ; j < size ; ++j )
				{
					if ( GetGeneType( seqs[ allOverlaps[j].seqIdx ].name ) != i )
						continue ;
					int l ;
					int seqIdx2 = allOverlaps[j].seqIdx ;
					if ( seqIdx2 == seqIdx ||
						!IsBetterGeneMatch( allOverlaps[j], geneOverlap[i], 0.95 )
						/*|| ( allOverlaps[j].readStart < geneOverlap[i].readStart - 20 
							|| allOverlaps[j].readStart > geneOverlap[i].readStart +  20 ) 
						|| ( allOverlaps[j].readEnd < geneOverlap[i].readEnd - 20 
							|| allOverlaps[j].readEnd > geneOverlap[i].readEnd +  20 )  */
					   )
						continue ;
					
					if ( IsSameGeneAllele( seqs[ seqIdx ].name, seqs[ seqIdx2 ].name ) )
						continue ;
					for ( l = 0 ; l < usedSeqIdx.Size() ; ++l )
					{
						if ( IsSameGeneAllele( seqs[ usedSeqIdx[l] ].name, seqs[ seqIdx2 ].name ) )
							break ;
					}
					if ( l < usedSeqIdx.Size() )
						continue ;

					++reportCnt ;
					sprintf( buffer + strlen( buffer ), ",%s(%d):(%d-%d):(%d-%d):%.2lf",
						seqs[ seqIdx2 ].name, seqs[ seqIdx2 ].consensusLen,
						allOverlaps[j].readStart, allOverlaps[j].readEnd, 
						allOverlaps[j].seqStart, allOverlaps[j].seqEnd, allOverlaps[j].similarity * 100 ) ;
					usedSeqIdx.PushBack( allOverlaps[j].seqIdx ) ;
					if ( reportCnt >= 2 )
						break ;
				}
			}
			else
			{
				sprintf( buffer + strlen( buffer ), " *" ) ;
			}
		}
		
		if ( cdr1 == NULL)
			sprintf( buffer + strlen( buffer), " CDR1(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR1(%d-%d):%.2lf=%s", cdr[0].readStart, cdr[0].readEnd, 
				cdr[0].similarity * 100, cdr1 ) ;
		
		if ( cdr2 == NULL)
			sprintf( buffer + strlen( buffer), " CDR2(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR2(%d-%d):%.2lf=%s", cdr[1].readStart, cdr[1].readEnd, 
				cdr[1].similarity * 100, cdr2 ) ;

		if ( cdr3 == NULL)
			sprintf( buffer + strlen( buffer), " CDR3(0-0):0.00=null" ) ;
		else
			sprintf( buffer + strlen( buffer), " CDR3(%d-%d):%.2lf=%s", cdr[2].readStart, cdr[2].readEnd, cdr3Score, cdr3 ) ;
		
		if ( cdr1 != NULL )
			delete[] cdr1 ;
		if ( cdr2 != NULL )
			delete[] cdr2 ;
		if ( cdr3 != NULL )
			delete[] cdr3 ;
		return 1 ;
	}
	
	// Use the refSet to annotate current set.
	void Annotate( SeqSet &refSet )
	{
		int i ;
		char *buffer = new char[10240] ;
		int seqCnt = seqs.size() ;
		struct _overlap geneOverlap[4];
		struct _overlap cdr[3] ;
		for ( i = 0 ; i < seqCnt  ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;
		
			free( seqs[i].name ) ;
			refSet.AnnotateRead( seqs[i].consensus, 2, geneOverlap, cdr, buffer ) ;
			seqs[i].name = strdup( buffer ) ;
		}

		delete[] buffer ;
	}
	
	
	void BreakFalseAssembly( std::vector<struct _Read> reads )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;
		std::vector<struct _overlap> *adj = new std::vector<struct _overlap>[ seqCnt ] ;
		//BuildBranchGraph( adj, 31 ) ;
		
		int readCnt = reads.size() ;
		
		// Mark whether to trust the branch overlap portion.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			std::vector<struct _overlap> overlaps ;
			int overlapCnt = GetOverlapsFromRead( reads[i].seq, 1, overlaps ) ;
			if ( overlapCnt == 0 )
				continue ;

			std::sort( overlaps.begin(), overlaps.end() ) ;
			
			struct _overlap extendedOverlap ;
			int len = strlen( reads[i].seq ) ;
			char *rc = strdup( reads[i].seq ) ;
			ReverseComplement( rc, reads[i].seq, len ) ;
			char *r = reads[i].seq ;
			if ( overlaps[0].strand == -1 )
				r = rc ;
			char *align = new char[ 2 * len + 2 ] ;	
			for ( j = 0 ; j < overlapCnt ; ++j )
			{
				if ( ExtendOverlap( r, len, seqs[ overlaps[j].seqIdx ], align, 
					overlaps[j], extendedOverlap ) == 1 )
				{
					int seqIdx = extendedOverlap.seqIdx ;
					int adjCnt = adj[ seqIdx ].size() ;
					for ( k = 0 ; k < adjCnt ; ++k )
					{
						if ( extendedOverlap.seqStart < adj[ seqIdx ][k].readStart &&
							extendedOverlap.seqEnd > adj[ seqIdx ][k].readEnd )
						{
							adj[ seqIdx ][k].similarity = 0 ;
						}
					}
					break ;
				}
			}
			free( rc ) ;
			delete[] align ;
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int size = adj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( adj[i][j].similarity < repeatSimilarity )
					adj[i][j].similarity = 0 ;

				//printf( "%d branch %d %d: %d %d %lf\n", i, j, adj[i][j].seqIdx, 
				//	adj[i][j].readStart, adj[i][j].readEnd, adj[i][j].similarity ) ;
			}
		}

		// Break up the Seq 
		int newSeqId = seqCnt ;
		KmerCode kmerCode( kmerLength ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int adjCnt = adj[i].size() ;
			// Find the regions.
			SimpleVector<struct _pair> breakCoords ;
			breakCoords.Reserve( adjCnt ) ;
			for ( k = 0 ; k < adjCnt ; ++k )
			{
				if ( adj[i][k].similarity <= 0 )
					continue ;
				
				struct _pair nb ;
				nb.a = adj[i][k].readStart ;
				nb.b = adj[i][k].readEnd ;
				if ( nb.a == 0 || nb.b == seqs[i].consensusLen - 1 )
					continue ;

				breakCoords.PushBack( nb ) ;
			}
			std::sort( breakCoords.BeginAddress(), breakCoords.EndAddress(), CompSortPairAInc ) ;	
			int size = breakCoords.Size() ;

			if ( size == 0 )
				continue ;
			k = 1 ;
			for ( j = 1 ; j < size ; ++j )
			{
				if ( breakCoords[j].a <= breakCoords[k - 1].b )
				{
					if ( breakCoords[j].b > breakCoords[k - 1].b )
						breakCoords[k - 1].b = breakCoords[j].b ;
				}
				else
				{
					breakCoords[k] = breakCoords[j] ;
					++k ;
				}
			}
			breakCoords.Resize( k ) ;
			
			int start, end ;
			for ( j = 0 ; j < k ; ++j )
				printf( "%d %s breakCoords %d: %d %d\n", i, seqs[i].name, j, breakCoords[j].a, breakCoords[j].b ) ;
			for ( j = 0 ; j <= k ; ++j )
			{
				if ( j == 0 )
					start = 0 ;
				else
					start = breakCoords[j - 1].a ;

				if ( j == k )
					end = seqs[i].consensusLen - 1 ;
				else
					end = breakCoords[j].b ;

				struct _seqWrapper ns ;
				ns.consensus = (char *)malloc( sizeof( char ) * ( end - start + 2 ) ) ;
				memcpy( ns.consensus, seqs[i].consensus + start, end - start + 1 ) ;
				ns.consensus[ end - start + 1 ] = '\0' ;
				ns.consensusLen = end - start + 1 ;
				ns.name = strdup( seqs[i].name ) ;
				ns.isRef = false ;
				
				ns.posWeight.Reserve( end - start + 1 ) ;
				int l ;
				for ( l = start ; l <= end ; ++l )
					ns.posWeight.PushBack( seqs[i].posWeight[l] ) ;	

				seqs.push_back( ns ) ;
				printf( "Break %d to %d.\n", i, seqs.size() - 1 ) ;
			}
			// Clean up original seq
			seqIndex.RemoveIndexFromRead( kmerCode, seqs[i].consensus, seqs[i].consensusLen, i, 0 ) ;
			free( seqs[i].consensus ) ;
			free( seqs[i].name ) ;
			seqs[i].consensus = seqs[i].name = NULL ;
			seqs[i].posWeight.Release() ;
		}
		//Clean( true ) ;	
	}

	// The readId is for the first mate, the other mate readId is readId+1.
	//   The main program should determine which one to use.
	void UpdateMateAdjGraph( int from, int fromStart, int fromEnd, int to, int toStart, int toEnd, 
			int readId, std::vector<struct _overlap> *mateAdj )
	{
		int j ;
		int size = mateAdj[from].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			if ( mateAdj[from][j].seqIdx == to )
			{
				if ( fromStart < mateAdj[from][j].readStart )
					mateAdj[from][j].readStart = fromStart ;
				if ( fromEnd > mateAdj[from][j].readEnd )
					mateAdj[from][j].readEnd = fromEnd ;
				if ( toStart < mateAdj[from][j].seqStart )
					mateAdj[from][j].seqStart = toStart ;
				if ( toEnd > mateAdj[from][j].seqEnd )
					mateAdj[from][j].seqEnd = toEnd ;
				++mateAdj[from][j].matchCnt ;
				mateAdj[from][j].info->PushBack( readId ) ;
				break ;
			}
		}
		if ( j >= size )
		{
			struct _overlap na ;
			na.seqIdx = to ;
			na.readStart = fromStart ;
			na.readEnd = fromEnd ;
			na.seqStart = toStart ;
			na.seqEnd = toEnd ;
			na.matchCnt = 1 ;
			na.info = new SimpleVector<int> ;
			na.info->PushBack( readId ) ;
			mateAdj[from].push_back( na ) ;
		}
	}

	// For those seqs with single read support, use its mate pair's sequence to decide extension
	//    which can use a much smaller overlap criterion to add new terms to branch graph.
	void AugmentBranchGraphByMate( std::vector<struct _overlap> *branchAdj, std::vector<struct _overlap> *prevAdj, 
			std::vector<struct _overlap> *nextAdj, std::vector<struct _assignRead> &reads )
	{
		int i, j, k ;
		int seqCnt = seqs.size() ;

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			if ( prevAdj[i].size() + nextAdj[i].size() != 1 )
				continue ;
			struct _seqWrapper &seq = seqs[i] ;
			/*for ( j = 0 ; j < seq.consensusLen ; ++j )
			{
				if ( seq.posWeight[j].Sum() != 1 )
					break ;
			}
			if ( j < seq.consensusLen )
				continue ;*/		
			int bsize = branchAdj[i].size() ;
			int readIdx, mateReadIdx ;
			char *fr ; // first read
			char *sr ; // second read;
			int flen, slen ; // their lengths

			// For overlaps
			int offset = -1 ;
			int overlapSize = -1 ;
			int bestMatchCnt = -1 ;

			if ( prevAdj[i].size() == 1 )
			{
				// Extend left
				for ( j = 0 ; j < bsize ; ++j )
				{
					if ( branchAdj[i][j].seqIdx == prevAdj[i][0].seqIdx )
						break ;
				}

				if ( j < bsize )
					continue ;

				//if ( prevAdj[i][0].info->Size() == 0 )
				//	continue ;
				int size = prevAdj[i][0].info->Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					readIdx = prevAdj[i][0].info->Get(j) ;
					mateReadIdx = readIdx + 1 ;
					if ( reads[ readIdx ].overlap.seqIdx != i )
					{
						int tmp = readIdx ;
						readIdx = mateReadIdx ;
						mateReadIdx = tmp ;
					}

					flen = strlen( reads[ mateReadIdx ].read  ) ;
					slen = strlen( reads[ readIdx ].read ) ;
					//if ( slen != seq.consensusLen || flen != seqs[ prevAdj[i][0].seqIdx ].consensusLen  )
					//	continue ;
					if ( reads[ readIdx ].overlap.seqStart > 3 )
						continue ;
					
					fr = strdup( reads[ mateReadIdx ].read ) ;
					sr = strdup( reads[ readIdx ].read ) ;
					int minOverlap = ( flen + slen) / 20 ;
					if ( minOverlap >= 31 )
						minOverlap = 31 ;

					if ( reads[ mateReadIdx ].overlap.strand == -1 )
						ReverseComplement( fr, reads[ mateReadIdx ].read, flen ) ;
					if ( reads[ readIdx ].overlap.strand == -1 )
						ReverseComplement( sr, reads[ readIdx ].read, slen ) ;
					overlapSize = AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, offset, bestMatchCnt ) ;
				
					if ( overlapSize == -1 )
					{
						free( fr ) ; free( sr ) ;
					}
					else
						break ;
				}
			}
			else
			{
				// Extend right
				for ( j = 0 ; j < bsize ; ++j )
				{
					if ( branchAdj[i][j].seqIdx == nextAdj[i][0].seqIdx )
						break ;
				}
				if ( j < bsize )
					continue ;

				//if ( nextAdj[i][0].info->Size() == 0 )
				//	continue ;
				int size = nextAdj[i][0].info->Size() ;
				for ( j = 0 ; j < size ; ++j )
				{
					readIdx = nextAdj[i][0].info->Get(j) ;
					mateReadIdx = readIdx + 1 ;

					if ( reads[ readIdx ].overlap.seqIdx != i )
					{
						int tmp = readIdx ;
						readIdx = mateReadIdx ;
						mateReadIdx = tmp ;
					}

					flen = strlen( reads[ readIdx ].read ) ;
					slen = strlen( reads[ mateReadIdx ].read ) ;
					//if ( flen != seq.consensusLen /*|| slen != seqs[ nextAdj[i][0].seqIdx ].consensusLen */)
					if ( reads[ readIdx ].overlap.seqEnd < seq.consensusLen - 4 )
						continue ;
					int minOverlap = ( flen + slen) / 20 ;
					if ( minOverlap >= 31 )
						minOverlap = 31 ;

					fr = strdup( reads[ readIdx ].read ) ;
					sr = strdup( reads[ mateReadIdx ].read ) ;

					if ( reads[ readIdx ].overlap.strand == -1 )
						ReverseComplement( fr, reads[ readIdx ].read, flen ) ;
					if ( reads[ mateReadIdx ].overlap.strand == -1 )
						ReverseComplement( sr, reads[ mateReadIdx ].read, slen ) ;
					overlapSize = AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, offset, bestMatchCnt ) ;
					
					if ( overlapSize == -1 )
					{
						free( fr ) ; free( sr ) ;
					}
					else
						break ;
				}
			}
			

			if ( overlapSize == -1 ) // Ambiguous overlap or no overlap
			{
				//free( fr ) ; free( sr ) ;
				continue ;
			}
			bestMatchCnt *= 2 ;
			if ( prevAdj[i].size() == 1 )
			{
				struct _overlap no ;
				no.seqIdx = prevAdj[i][0].seqIdx ;
				no.readStart = reads[ readIdx ].overlap.seqStart ;
				no.readEnd = no.readStart + overlapSize - 1 ;
				no.seqStart = reads[ mateReadIdx ].overlap.seqStart + offset ;
				no.seqEnd = no.seqStart + overlapSize - 1 ;
				no.matchCnt = bestMatchCnt ;
				no.similarity = bestMatchCnt / (2.0 * overlapSize ) ;
				branchAdj[i].push_back( no ) ;
				
				//printf( "<=%d %d\n%s\n%s\n\n", no.seqIdx, i, fr, sr ) ;
			}
			else
			{
				struct _overlap no ;
				no.seqIdx = nextAdj[i][0].seqIdx ;
				no.readStart = reads[ readIdx ].overlap.seqStart + offset ;
				no.readEnd = no.readStart + overlapSize - 1 ;
				no.seqStart = reads[ mateReadIdx ].overlap.seqStart ;
				no.seqEnd = no.seqStart + overlapSize - 1 ;
				no.matchCnt = bestMatchCnt ;
				no.similarity = bestMatchCnt / (2.0 * overlapSize) ;
				branchAdj[i].push_back( no ) ;
				//printf( ">=%d %d\n%s\n%s\n\n", i, no.seqIdx, fr, sr ) ;
			}
			free( fr ) ; free( sr ) ;
		}
	}
	
	// return: 0: no extension. 1: end-to-inside extension. 2: end-to-end extension
	int GetExtendSeqCoord( int from, struct _overlap mateInfo, int direction, 
			std::vector<struct _overlap> *branchAdj, bool aggressive, struct _overlap &coord )
	{
		int i, k ;
		int adjCnt = branchAdj[from].size() ;	
		int to = mateInfo.seqIdx ;

		coord.seqIdx = -1 ;
		int overhangSize = 5 ; // Allowing the end to have this much random extension in the raw assembly.
		
		for ( i = 0 ; i < adjCnt ; ++i )
		{
			if ( direction == 1 )
			{
				if ( branchAdj[from][i].seqIdx == to 
						&& branchAdj[from][i].readEnd >= seqs[from].consensusLen - overhangSize ) 
						//&& branchAdj[from][i].seqStart <= 5 - 1 )
					break ;
			}
			else if ( direction == -1 )
			{
				if ( branchAdj[from][i].seqIdx == to 
						&& branchAdj[from][i].readStart <= overhangSize - 1 ) 
						//&& branchAdj[from][i].seqEnd <= seqs[to].consensusLen - 5 )
					break ;
			}
		}


		if ( i >= adjCnt )
			return 0 ;
		k = i ;
		
		if ( direction == 1 && mateInfo.seqEnd <= branchAdj[from][k].seqEnd )
			return 0 ;
		else if ( direction == -1 && mateInfo.seqStart >= branchAdj[from][k].seqStart )
			return 0 ;

		/*for ( i = 0 ; i < adjCnt ; ++i )
		{
		}*/
		
		coord.seqIdx = to ;
		coord.matchCnt = branchAdj[from][k].readEnd - branchAdj[from][k].readStart + 1 ; // Record the length of the branch overlap.
		int ret = 1 ;
		
		if ( direction == 1 )
		{
			coord.readStart = 0 ;
			coord.readEnd = branchAdj[from][k].readEnd ;
		}
		else
		{
			coord.readStart = branchAdj[from][k].readStart ;
			coord.readEnd = seqs[from].consensusLen - 1 ;
		}

		if ( direction == 1 )
		{
			coord.seqStart = branchAdj[from][k].seqEnd + 1 ;
			if ( aggressive )
				coord.seqEnd = seqs[to].consensusLen - 1 ;
			else
				coord.seqEnd = mateInfo.seqEnd ;

			if ( branchAdj[from][k].seqStart <= overhangSize - 1 )
				ret = 2 ;
		}
		else
		{
			if ( aggressive )
				coord.seqStart = 0 ;
			else
			{
				//printf( "%d %d: %d %d\n", from, to, mateInfo.seqStart, branchAdj[from][k].seqStart ) ;
				coord.seqStart = mateInfo.seqStart ;
			}
			coord.seqEnd = branchAdj[from][k].seqStart - 1 ;

			if ( branchAdj[from][k].seqEnd >= seqs[ branchAdj[from][k].seqIdx ].consensusLen - overhangSize )
				ret = 2 ;
		}

		return ret ;
	}
	
	// Extend seq where one mate is aligned the other mate is missing, and the mates overlap with each other.
	void ExtendSeqFromMissingOverlapMate( std::vector< struct _assignRead> reads )
	{
		int i, j, k ;
		int readCnt = reads.size() ;
		int seqCnt = seqs.size() ;
		for ( i = 0 ; i < readCnt ; ++i )
		{
			bool paired = false ;
			if ( i < readCnt - 1 && !strcmp( reads[i].id, reads[i + 1].id ) )
				paired = true ;
			if ( !paired )
				continue ;

			if ( ( reads[i].overlap.seqIdx != -1 && reads[i + 1].overlap.seqIdx != -1 )
				|| ( reads[i].overlap.seqIdx == -1 && reads[i + 1].overlap.seqIdx == -1 ) )
			{
				++i ;
				continue ;
			}
			
			int anchorId = i ;
			int anchorSeqIdx = reads[i].overlap.seqIdx ;
			int hangId = i + 1 ;
			int asLen, hsLen ;
			if ( reads[i].overlap.seqIdx == -1 )		
			{
				anchorId = i + 1 ;
				anchorSeqIdx = reads[i + 1].overlap.seqIdx ;
				hangId = i ;
			}
		
			struct _seqWrapper &seq = seqs[ anchorSeqIdx ] ;
			//printf( "%d %s\n%s\n", anchorSeqIdx, seq.consensus, reads[ anchorId ].read ) ;
			if ( ( reads[ anchorId ].overlap.strand == 1 && strcmp( seq.consensus, reads[ anchorId ].read ) ) 
				|| ( reads[ anchorId ].overlap.strand == -1 && !IsReverseComplement( seq.consensus, reads[ anchorId ].read ) ) )
			{
				++i ;
				continue ;
			}
			asLen = strlen( reads[ anchorId ].read ) ;
			hsLen = strlen( reads[ hangId ].read ) ;
			
			char *r = strdup( reads[ hangId ].read ) ;
			char *fr, *sr ;
			int direction = 0 ;
			int flen, slen ;
			if ( reads[ anchorId ].overlap.strand == 1 )
			{
				// Extend towards right
				ReverseComplement( r, reads[ hangId ].read, hsLen ) ;

				fr = reads[ anchorId ].read ;
				flen = asLen ;
				sr = r ;
				slen = hsLen ;
				direction = 1 ;
			}
			else
			{
				fr = r ;
				flen = hsLen ;
				sr = seq.consensus ;
				slen = asLen ;
				direction = -1 ;
			}
			
			int offset ;
			int matchCnt ;
			int minOverlap = ( flen + slen) / 20 ;
			if ( minOverlap >= 31 )
				minOverlap = 31 ;
			//printf( "overlap test %d %d %d\n%s\n%s\n", flen, slen, minOverlap, fr, sr ) ;
			if ( AlignAlgo::IsMateOverlap( fr, flen, sr, slen, minOverlap, 
				offset, matchCnt ) == -1 )
			{
				++i ;
				free( r ) ;
				continue ;
			}

			if ( flen - offset >= slen ) // sr is contained in fr
			{
				++i ;
				free( r ) ;
				continue ;
			}

			int newConsensusLen = offset + slen ; 
			char *newConsensus = (char*)malloc( sizeof( char ) * ( newConsensusLen + 1 ) ) ;
			memcpy( newConsensus, fr, offset ) ;
			memcpy( newConsensus + offset, sr, slen ) ;
			newConsensus[newConsensusLen] = '\0' ;
			
			// Update the pos weight.
			if ( direction == 1 )
			{
				// Append
				seq.posWeight.ExpandTo( newConsensusLen ) ;
				seq.posWeight.SetZero( flen, newConsensusLen - flen ) ;
				UpdatePosWeightFromRead( seq.posWeight, offset, sr ) ;
				/*for ( j = 0 ; j < slen ; ++j )
				{
					if ( sr[j] == 'N' )
						continue ;
					++seq.posWeight[j + offset].count[ nucToNum[ sr[j] - 'A' ] ] ;
				}*/
			}
			else
			{
				seq.posWeight.ShiftRight( offset ) ;
				seq.posWeight.SetZero( 0, offset ) ;
				UpdatePosWeightFromRead( seq.posWeight, 0, fr ) ;
				/*for ( j = 0 ; j < flen ; ++j )
				{
					if ( fr[j] == 'N' )
						continue ;
					++seq.posWeight[j].count[ nucToNum[ fr[j] - 'A' ] ] ;
				}*/
			}
			
			free( r ) ;
			//printf( "pass\n%s\n%s\n%d %d\n", seq.consensus, newConsensus, seq.consensusLen, newConsensusLen ) ;
			//fflush( stdout ) ;
			free( seq.consensus ) ;
			seq.consensus = newConsensus ;
			seq.consensusLen = newConsensusLen ;
			seq.minLeftExtAnchor = seq.minRightExtAnchor = 0 ;
			++i ;
		}
	}

	// Use this set of reads to extend,rearrange the seq 
	void ExtendSeqFromReads( std::vector<struct _assignRead> &reads, int leastOverlapLen )
	{
		int i, j, k ;
		int readCnt = 0 ;
		int seqCnt = seqs.size() ;
		double backupNovelSeqSimilarity = novelSeqSimilarity ;

		novelSeqSimilarity = 1.00 ;
		int ret = 0 ;

		std::vector<struct _overlap> *branchAdj = new std::vector<struct _overlap>[ seqCnt ] ;

		// Mate-pair support graph. In here, start,end represnts the rough range with read support
		//	and matchCnt represent how many mates support this connection.
		std::vector<struct _overlap> *nextAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		std::vector<struct _overlap> *prevAdj = new std::vector<struct _overlap>[ seqCnt ] ; 
		readCnt = reads.size() ;
		
		SimpleVector<bool> useInBranch ;
		useInBranch.ExpandTo( seqCnt ) ;
		useInBranch.SetZero( 0, seqCnt ) ;

		// Directly use overlapped mate pairs for extension on 
		//   singleton assembly where the other mate has no perfect alignment.
		// Since these reads will not be applied on more sophisticated extension, it is fine 
		//   to make it an independent component.
		std::sort( reads.begin(), reads.end(), CompSortAssignedReadById ) ;
		//ExtendSeqFromMissingOverlapMate( reads ) ;

		// Then do more sophisticated extension
		// Build the mate adj graph.
		for ( i = 0 ; i < readCnt ; ++i )
		{
			bool paired = false ;
			if ( i < readCnt - 1 && !strcmp( reads[i].id, reads[i + 1].id ) )
				paired = true ; 
			
			if ( paired )
			{
				if ( reads[i].overlap.seqIdx == -1 || reads[i + 1].overlap.seqIdx == -1 
					|| reads[i].overlap.strand == reads[i + 1].overlap.strand )				
				{
					++i ;
					continue ;
				}
				int from, to ;
				int fromStart, fromEnd, toStart, toEnd ;
				struct _overlap extendedOverlap = reads[i].overlap ;
				struct _overlap mateExtendedOverlap = reads[i + 1].overlap ;
								
				if ( extendedOverlap.strand == 1 )
				{
					from = extendedOverlap.seqIdx ;
					fromStart = extendedOverlap.seqStart ;
					fromEnd = extendedOverlap.seqEnd ;

					to = mateExtendedOverlap.seqIdx ;
					toStart = mateExtendedOverlap.seqStart ;
					toEnd = mateExtendedOverlap.seqEnd ;
				}
				else
				{
					to = extendedOverlap.seqIdx ;
					toStart = extendedOverlap.seqStart ;
					toEnd = extendedOverlap.seqEnd ;

					from = mateExtendedOverlap.seqIdx ;
					fromStart = mateExtendedOverlap.seqStart ;
					fromEnd = mateExtendedOverlap.seqEnd ;
				}
			
				useInBranch[from] = true ;
				useInBranch[to] = true ;
				UpdateMateAdjGraph( from, fromStart, fromEnd, to, toStart, toEnd, i, nextAdj ) ;
				UpdateMateAdjGraph( to, toStart, toEnd, from, fromStart, fromEnd, i, prevAdj ) ;
				
				//printf( "%d(%d %d) %d(%d %d)\n", 
				//	from, fromStart, fromEnd,
				//	to, toStart, toEnd ) ;
			}
			else
			{
				// Single-end reads, not sure how to use them for now.
				// 	may be look for fragmented secondary alignment?
			}
			
			if ( paired )
				++i ;
		
		}

		BuildBranchGraph( branchAdj, leastOverlapLen, useInBranch ) ;
		//AugmentBranchGraphByMate( branchAdj, prevAdj, nextAdj, reads ) ;
		
		// Each Seq will be extend to left and right one step.
		SimpleVector<int> toRemoveSeqIdx ;
		toRemoveSeqIdx.Reserve( seqCnt ) ;
		
		// Find which seq id extend to.
		SimpleVector<struct _pair> matePrevNext ; 
		matePrevNext.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevAdjCnt = prevAdj[i].size() ;
			int prevTag = -1 ;
			int max = -1 ;
			for ( j = 0 ; j < prevAdjCnt ; ++j )
			{
				if ( prevAdj[i][j].seqIdx == i )
					continue ;

				if ( prevAdj[i][j].matchCnt > max )
				{
					prevTag = j ;
					max = prevAdj[i][j].matchCnt ;
				}
				else if ( prevAdj[i][j].matchCnt == max )
				{
					if ( prevAdj[i][j].seqEnd - prevAdj[i][j].seqStart > 
							prevAdj[i][ prevTag ].seqEnd - prevAdj[i][prevTag].seqStart )
						prevTag = j ;
				}
#ifdef DEBUG				
				printf( "<= %d: %d %d. %d %d %d %d\n", i, prevAdj[i][j].seqIdx, prevAdj[i][j].matchCnt, 
						prevAdj[i][j].readStart, prevAdj[i][j].readEnd,
						prevAdj[i][j].seqStart, prevAdj[i][j].seqEnd ) ;
#endif
			}

			int nextAdjCnt = nextAdj[i].size() ;
			int nextTag = -1 ;
			max = -1 ;
			for ( j = 0 ; j < nextAdjCnt ; ++j )
			{
				if ( nextAdj[i][j].seqIdx == i )
					continue ;

				if ( nextAdj[i][j].matchCnt > max )
				{
					nextTag = j ;
					max = nextAdj[i][j].matchCnt ;
				}
				else if ( nextAdj[i][j].matchCnt == max )
				{
					if ( nextAdj[i][j].seqEnd - nextAdj[i][j].seqStart > 
							nextAdj[i][ nextTag ].seqEnd - nextAdj[i][nextTag].seqStart )
						nextTag = j ;
				}

#ifdef DEBUG
				printf( "=> %d: %d %d. %d %d %d %d\n", i, nextAdj[i][j].seqIdx, nextAdj[i][j].matchCnt, 
						nextAdj[i][j].readStart, nextAdj[i][j].readEnd,
						nextAdj[i][j].seqStart, nextAdj[i][j].seqEnd ) ;
#endif
			}
			matePrevNext[i].a = prevTag ;
			matePrevNext[i].b = nextTag ;
		}

		// Figure out the seqs whose extension will create repeat sequence
		// 	e.g.: end-to-end extension from two anchor seq.
		SimpleVector<struct _pair> extensionType ;
		SimpleVector<int> uniqueSuccessorOf ;
		extensionType.ExpandTo( seqCnt ) ;
		uniqueSuccessorOf.ExpandTo( seqCnt ) ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;

			struct _overlap leftExtend, rightExtend ;
			leftExtend.seqIdx = -1 ;
			rightExtend.seqIdx = -1 ;
			extensionType[i].a = extensionType[i].b = 0 ;
			if ( prevTag >= 0 )
			{
				extensionType[i].a = GetExtendSeqCoord( i, prevAdj[i][ prevTag ], -1, branchAdj, false, leftExtend ) ;
				if ( leftExtend.seqIdx == -1 )
					matePrevNext[i].a = -1 ;

			}
			if ( nextTag >= 0 )
			{
				extensionType[i].b = GetExtendSeqCoord( i, nextAdj[i][ nextTag ], 1, branchAdj, false, rightExtend ) ;
				if ( rightExtend.seqIdx == -1 )
					matePrevNext[i].b = -1 ;
			}
		}

		// Rescue some partial extension. e.g.: 
		// A |;  A<=B. A could not extend right due to no overlap, then we will put A=>B back.
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;

			if ( prevTag >= 0 )
			{
				int seqIdx = prevAdj[i][ prevTag ].seqIdx ;
				if ( matePrevNext[ seqIdx ].b == -1 && extensionType[i].a == 2 )
				{
					// rescue the connection.
					int size = nextAdj[ seqIdx ].size() ;
					for ( j = 0 ; j < size ; ++j )
						if ( nextAdj[ seqIdx ][j].seqIdx == i )
						{
							struct _overlap tmp ;
							extensionType[ seqIdx ].b = GetExtendSeqCoord( seqIdx, nextAdj[ seqIdx ][j], 
										1, branchAdj, false, tmp ) ;
							if ( extensionType[ seqIdx ].b == 2 )
								matePrevNext[ seqIdx ].b = j ;
							break ;
						}
				}
			}

			if ( nextTag >= 0 )
			{
				int seqIdx = nextAdj[i][ nextTag ].seqIdx ;
				if ( matePrevNext[ seqIdx ].a == -1 && extensionType[i].b == 2 )
				{
					int size = prevAdj[seqIdx].size() ;
					for ( j = 0 ; j < size; ++j )
						if ( prevAdj[ seqIdx ][j].seqIdx == i )
						{
							struct _overlap tmp ;
							extensionType[ seqIdx ].a = GetExtendSeqCoord( seqIdx, prevAdj[ seqIdx ][j], 
										-1, branchAdj, false, tmp ) ;
							if ( extensionType[ seqIdx ].a == 2 )
								matePrevNext[ seqIdx ].a = j ;
							break ;
						}
				}
			}
		}
		
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			uniqueSuccessorOf[i] = -1 ;
			
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ;
			//if ( i == 47 )
			//	fprintf( stderr, "hi %d %d %d %d\n", prevTag, nextTag,
			//		prevAdj[i][prevTag].seqIdx, nextAdj[i][nextTag].seqIdx ) ;
			//if ( prevTag >= 0 && nextTag >= 0 )
			//	continue ;
			/*if ( nextTag >= 0 )
			{
				// Keep the second part.
				if ( extensionType[i].b == 2 )
				{
					int seqIdx = nextAdj[i][nextTag].seqIdx ;
					if ( matePrevNext[ seqIdx ].a >= 0 
						&& prevAdj[ seqIdx ][ matePrevNext[ seqIdx ].a ].seqIdx == i 
						&& extensionType[ seqIdx ].a == 2 )
					{
						directlyFilter[i] = true ;
					}
				}
			}*/

			if ( prevTag >= 0 )
			{
				if ( extensionType[i].a == 2 )
				{
					int seqIdx = prevAdj[i][ prevTag ].seqIdx ;
					if ( matePrevNext[seqIdx].b >= 0 
						&& nextAdj[ seqIdx ][ matePrevNext[ seqIdx ].b ].seqIdx == i 
						&& extensionType[ seqIdx ].b == 2 )
					{
						uniqueSuccessorOf[i] = seqIdx ;
					}
				}
			}
		}
		

		// Filter some middle part, if its two extension are connected.
		/*for ( i = 0 ; i < seqCnt ; ++i )
		{
			int prevTag = matePrevNext[i].a ;
			int nextTag = matePrevNext[i].b ; 
			if ( directlyFilter[i] || prevTag < 0 || nextTag < 0 || extensionType[i].a != 2 
				|| extensionType[i].b != 2 )
				continue ;
			
			int prevSeqIdx = prevAdj[i][ prevTag ].seqIdx ;
			int nextSeqIdx = nextAdj[i][ nextTag ].seqIdx ;
			if ( ( matePrevNext[ prevSeqIdx ].b >= 0 
					&& nextAdj[ prevSeqIdx ][ matePrevNext[ prevSeqIdx ].b ].seqIdx == nextSeqIdx ) 
				|| ( matePrevNext[ nextSeqIdx ].a >= 0 
					&& prevAdj[ nextSeqIdx ][ matePrevNext[ nextSeqIdx ].a ].seqIdx == prevSeqIdx ) )
				directlyFilter[i] = true ;
		}*/

		// Do the extesion
		SimpleVector<int> chain ;
		SimpleVector<int> offset ;
		SimpleVector<struct _pair> range ; // the range of a seq that needs to be copied in the new seq.
		SimpleVector<int> origRangeB ;
		SimpleVector<struct _pair> shiftSeq ; // where the seq moved to after extesion. a-new seqIdx, b-offset in new seq.
					              // Note that the ending of a seq might be trimmed when extension,
						      //  the b(offset) will regard those trimmed part existing, and hence
						      //  starts a bit earlier than the portion that is actually used.
		chain.Reserve( seqCnt ) ;	    
		offset.Reserve( seqCnt ) ;
		range.Reserve( seqCnt ) ;
		origRangeB.Reserve( seqCnt ) ; // the range.b that are not adjusted.

		shiftSeq.ExpandTo( seqCnt ) ; 
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			shiftSeq[i].a = i ;
			shiftSeq[i].b = 0 ;
		}

		for ( i = 0 ; i < seqCnt ; ++i )
		{
			//if ( inCnt[i] > 0 )
			//	continue ;
			// Each seq will try to extend to each direction once.
			if ( uniqueSuccessorOf[i] != -1 )
			{
				toRemoveSeqIdx.PushBack( i ) ;
				continue ;
			}
			
			int first = i ;
			int firstPrevTag = matePrevNext[i].a ;
			int last = i ;
			int lastNextTag = matePrevNext[i].b ;
			
			chain.Clear() ;
			chain.PushBack( i ) ;
			while ( 1 )
			{
				if ( lastNextTag >= 0 && uniqueSuccessorOf[ nextAdj[last][ lastNextTag ].seqIdx ] == last )
				{
					last = nextAdj[last][ lastNextTag ].seqIdx ;
					lastNextTag = matePrevNext[last].b ;
					chain.PushBack( last ) ;
				}
				else
					break ;
			}
			int chainSize = chain.Size() ;
			// Compute the length of new seq.
			int newConsensusLen = 0 ;
			offset.Clear() ;
			range.Clear() ; 
			origRangeB.Clear() ;
			struct _overlap leftMostExtend, rightMostExtend ;
			leftMostExtend.seqIdx = -1 ;
			rightMostExtend.seqIdx = -1 ;
			offset.ExpandTo( chainSize ) ;
			range.ExpandTo( chainSize ) ;
			origRangeB.ExpandTo( chainSize ) ;
			for ( j = 0 ; j < chainSize ; ++j )
			{
				struct _overlap leftExtend, rightExtend ;
				int prevTag = matePrevNext[ chain[j] ].a ;
				int nextTag = matePrevNext[ chain[j] ].b ;
				leftExtend.seqIdx = -1 ;
				rightExtend.seqIdx = -1 ;
				
				if ( prevTag >= 0 )
					GetExtendSeqCoord( chain[j], prevAdj[ chain[j] ][ prevTag ], -1, branchAdj, j > 0, leftExtend ) ;
				if ( nextTag >= 0 )
					GetExtendSeqCoord( chain[j], nextAdj[ chain[j] ][ nextTag ], 1, branchAdj, j < chainSize - 1, rightExtend ) ;
				
				if ( j == 0 && leftExtend.seqIdx != -1 )
				{
					newConsensusLen += leftExtend.seqEnd - leftExtend.seqStart + 1 ;
					leftMostExtend = leftExtend ;
				}
				offset[j] = newConsensusLen ;
				if ( leftExtend.seqIdx != -1 )
					range[j].a = leftExtend.readStart ;
				else
					range[j].a = 0 ;

				if ( rightExtend.seqIdx != -1 )
					range[j].b = rightExtend.readEnd ;
				else
					range[j].b = seqs[ chain[j] ].consensusLen - 1 ;
				origRangeB[j] = range[j].b ;

				if ( j < chainSize - 1 )
				{
					// For the middle ones, we only keep the left overlap.
					/*struct _overlap tmp ;
					GetExtendSeqCoord( chain[j + 1], prevAdj[ chain[j + 1] ][ matePrevNext[ chain[j + 1] ].a ], 
						-1, branchAdj, tmp ) ;
					
					if ( tmp.seqIdx != -1 )
					{
						range[j].b = tmp.seqEnd ;
						if ( range[j].b < range[j].a )
							range[j].b = range[j].a - 1 ;
					}
					else
					{
						// should never get here.
						chainSize = j + 1 ;
					}*/
					range[j].b -= rightExtend.matchCnt ;
					if ( range[j].b < range[j].a )
						range[j].b = range[j].a - 1 ;
				} 

				newConsensusLen += range[j].b - range[j].a + 1 ;
				if ( j == chainSize - 1 && rightExtend.seqIdx != -1 )
				{
					newConsensusLen += rightExtend.seqEnd - rightExtend.seqStart + 1 ;
					rightMostExtend = rightExtend ;
				}
			}
			if ( newConsensusLen == seqs[i].consensusLen ) // no extension.
			{
				continue ;
			}
			//printf( "%d %d %d\n", chainSize, seqs[i].consensusLen, newConsensusLen ) ;	
			char *newConsensus = ( char * )malloc( sizeof( char ) * ( newConsensusLen + 1 ) ) ;
			// Put in the seqs.
			if ( leftMostExtend.seqIdx != -1 )
				memcpy( newConsensus, seqs[ leftMostExtend.seqIdx ].consensus + leftMostExtend.seqStart,
					leftMostExtend.seqEnd - leftMostExtend.seqStart + 1 ) ;
			for ( j = 0 ; j < chainSize ; ++j )
			{
				memcpy( newConsensus + offset[j], seqs[ chain[j] ].consensus + range[j].a,
					range[j].b - range[j].a + 1 ) ;
			}
			if ( rightMostExtend.seqIdx != -1 )
			{
				memcpy( newConsensus + offset[j - 1] + range[j - 1].b - range[j - 1].a + 1,
					seqs[ rightMostExtend.seqIdx ].consensus + rightMostExtend.seqStart,
					rightMostExtend.seqEnd - rightMostExtend.seqStart + 1 ) ;
			}
			newConsensus[newConsensusLen] = '\0' ;

			struct _seqWrapper ns ;
			ns.isRef = false ;
			ns.consensus = newConsensus ;
			ns.consensusLen = newConsensusLen ;
			
			ns.name = strdup( seqs[i].name ) ;
			ns.posWeight.ExpandTo( newConsensusLen ) ;
			ns.posWeight.SetZero( 0, newConsensusLen ) ;

			// Assign the posWeight of the core part.	
			for ( j = 0 ; j < chainSize ; ++j )
			{
				int l ;
				for ( l = range[j].a ; l <= origRangeB[j] && offset[j] + l - range[j].a < newConsensusLen ; ++l )
				{
					ns.posWeight[ offset[j] + l - range[j].a ] += seqs[ chain[j] ].posWeight[l] ;
				}
			}

			// Adjust the posWeight for the overhang part
			if ( leftMostExtend.seqIdx != -1 )
			{
				int from = leftMostExtend.seqIdx ;
				int to = chain[0] ;
				
				int size = nextAdj[from].size() ;
				for ( j = 0 ; j < size ; ++j )
					if ( nextAdj[from][j].seqIdx == to )
						break ;
				
				if ( j < size )
				{
					// Should always get here
					SimpleVector<int> &readIdx = *( nextAdj[from][j].info ) ; 
					size = readIdx.Size() ;
					int l ;
					for ( l = 0 ; l < size ; ++l )
					{
						int ridx ;
						if ( reads[ readIdx[l] ].overlap.seqIdx == from )
							ridx = readIdx[l] ;
						else
							ridx = readIdx[l] + 1 ;
						
						if ( reads[ ridx ].overlap.seqEnd > leftMostExtend.seqEnd + leftMostExtend.matchCnt )
							continue ;

						int m, rm ;
						for ( m = reads[ ridx ].overlap.seqStart, rm = 0 ; 
							m <= reads[ ridx ].overlap.seqEnd ; ++m, ++rm )
						{
							if ( reads[ ridx ].read[rm] != 'N' )
							{
								if ( m < newConsensusLen )
									++ns.posWeight[m - leftMostExtend.seqStart].
										count[ nucToNum[ reads[ ridx ].read[rm] - 'A' ] ] ;
								
								if ( shiftSeq[from].b + m >= 0 && 
									shiftSeq[from].b + m < seqs[ shiftSeq[from].a ].consensusLen )
									--seqs[ shiftSeq[from].a ].posWeight[ shiftSeq[from].b + m ].
												count[ nucToNum[ reads[ ridx ].read[rm] - 'A' ] ] ;
							}
						}
					}
				}
			}

			if ( rightMostExtend.seqIdx != -1 )
			{
				int from = chain[ chainSize - 1 ] ;
				int to = rightMostExtend.seqIdx ;
				
				int size = nextAdj[from].size() ;
				for ( j = 0 ; j < size ; ++j )
					if ( nextAdj[from][j].seqIdx == to )
						break ;
				
				if ( j < size )
				{
					// Should always get here
					SimpleVector<int> &readIdx = *( nextAdj[from][j].info ) ; 
					size = readIdx.Size() ;
					int l ;
					int lastOffset = offset[chainSize - 1] + range[chainSize - 1].b - range[chainSize - 1].a + 1 ; 
					for ( l = 0 ; l < size ; ++l )
					{
						int ridx ;
						if ( reads[ readIdx[l] ].overlap.seqIdx == from )
							ridx = readIdx[l] + 1 ;
						else
							ridx = readIdx[l] ;
						
						if ( reads[ ridx ].overlap.seqStart < rightMostExtend.seqStart - rightMostExtend.matchCnt )
							continue ;

						int m ;
						int rm ;
						char *s = strdup( reads[ ridx ].read ) ;
						if ( reads[ ridx ].overlap.strand == -1 )
						{
							// should always get here
							ReverseComplement( s, reads[ ridx ].read, strlen( reads[ ridx ].read ) ) ;
						}
						for ( m = reads[ ridx ].overlap.seqStart, rm = 0 ; 
							m <= reads[ ridx ].overlap.seqEnd ; ++m, ++rm )
						{
							if ( s[rm] != 'N' )
							{
								int adjustM = m - rightMostExtend.seqStart + lastOffset ;  
								//printf( "%d %d %d. %c %c\n", m, rm, adjustM, ns.consensus[adjustM], reads[ ridx ].read[rm] ) ;
								if ( adjustM >= 0 && adjustM < newConsensusLen )
									++ns.posWeight[adjustM].count[ nucToNum[ s[rm] - 'A' ] ] ;
								if ( shiftSeq[to].b + m >= 0 && 
									shiftSeq[to].b + m < seqs[ shiftSeq[to].a ].consensusLen )
									--seqs[ shiftSeq[to].a ].posWeight[ shiftSeq[to].b + m].
												count[ nucToNum[ s[rm] - 'A' ] ] ;
							}
						}
						free( s ) ;
					}
				}

			}
			
			// For other not updated region, just assign a number there.
			for ( j = 0 ; j < newConsensusLen ; ++j )
				if ( newConsensus[j] != 'N' && 
					ns.posWeight[j].Sum() == 0 )
					ns.posWeight[j].count[ nucToNum[ newConsensus[j] - 'A' ] ] = 1 ;

			// Update the shift information
			for ( j = 0 ; j < chainSize ; ++j )
			{
				shiftSeq[ chain[j] ].a = seqs.size() ;
				shiftSeq[ chain[j] ].b = offset[j] - range[j].a ;
			}

#ifdef DEBUG	
			if ( leftMostExtend.seqIdx != -1 )
				printf( "left 0: %d %s\n", leftMostExtend.seqIdx, seqs[ leftMostExtend.seqIdx ].consensus ) ;

			for ( j = 0 ; j < chainSize ; ++j )
				printf( "chain %d: %d %s\n", j + 1, chain[j], seqs[ chain[j] ].consensus ) ;
			if ( rightMostExtend.seqIdx != -1 )
				printf( "right %d: %d %s\n", j + 1, rightMostExtend.seqIdx, seqs[ rightMostExtend.seqIdx ].consensus ) ;
			printf( "%d new %s\n", i, newConsensus) ;
			fflush( stdout ) ;
#endif 

			ns.minLeftExtAnchor = ns.minRightExtAnchor = 0 ;
			seqs.push_back( ns ) ;
			toRemoveSeqIdx.PushBack( i ) ;
		}
		
		int size ;
		delete[] branchAdj ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			size = nextAdj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( nextAdj[i][j].info != NULL )
					delete nextAdj[i][j].info ;
			}
			
			size = prevAdj[i].size() ;
			for ( j = 0 ; j < size ; ++j )
			{
				if ( prevAdj[i][j].info != NULL )
					delete prevAdj[i][j].info ;
			}

		}
		delete[] nextAdj ;
		delete[] prevAdj ;
		
		size = toRemoveSeqIdx.Size() ;
		for ( i = 0 ; i < size ; ++i )
			ReleaseSeq( toRemoveSeqIdx[i] ) ;	
		Clean( true ) ;
		
		// Recompute the seq id the read is assigned to.

		// Recompute the posweight that becomes negative 
		// since the alignment might be changed when adding and after adding states
		seqCnt = seqs.size() ;
		for ( i = 0 ; i < seqCnt ; ++i )
		{
			for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
			{
				int sum = 0 ;
				for ( k = 0 ; k < 4 ; ++k )
				{
					if ( seqs[i].posWeight[j].count[k] < 0 )
						seqs[i].posWeight[j].count[k] = 0 ;
					sum += seqs[i].posWeight[j].count[k] ;
				}	
				if ( sum == 0 && seqs[i].consensus[j] != 'N' )
				{
					seqs[i].posWeight[j].count[ nucToNum[ seqs[i].consensus[j] - 'A' ] ] = 1 ;
				}
			}
		}

		novelSeqSimilarity = backupNovelSeqSimilarity ;
	}

	void Output( FILE *fp )
	{
		int i, j, k ;
		int size = seqs.size() ;
		for ( i = 0 ; i < size ; ++i )
		{
			if ( seqs[i].isRef || seqs[i].consensus == NULL )
				continue ;

			fprintf( fp, ">assemble%d %s\n%s\n", i, seqs[i].name, seqs[i].consensus ) ;
			
			for ( k = 0 ; k < 4 ; ++k )
			{
				for ( j = 0 ; j < seqs[i].consensusLen ; ++j )
					fprintf( fp, "%d ", seqs[i].posWeight[j].count[k] ) ;
				fprintf( fp, "\n" ) ;
			}
		}
	}

	char *GetSeqName( int seqIdx )
	{
		return seqs[ seqIdx ].name ;
	}

} ;


#endif
