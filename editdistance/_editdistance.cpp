// -------
// License
// -------
//
// It is released under the MIT license.
//
//     Copyright (c) 2013 Hiroyuki Tanaka
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//     The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <bitset>
#include <array>

#include "./_editdistance.h"

using namespace std;

enum class MAPOPTION {
	MAP,
	N1,
	N2
};


string to_str(uint64_t const &n) {
	string s;
	for(size_t i = 0; i < 64; ++i) {
		if(n & (1L << i)) s = "1" + s;
		else s = "0" + s;
	}
	return s;
}


template<typename T, typename V>
unsigned int edit_distance_bp(T &cmap, V const &vec, unsigned int m) {
	unsigned int D = m;
	uint64_t D0, HP, HN, VP = 0, VN = 0;
	uint64_t top = (1L << (m - 1));
	for(size_t i = 0; i < m; ++i) VP |= (1L << i);
	for(size_t i = 0; i < vec.size(); ++i) {
		uint64_t PM = cmap[vec[i]];
		D0 = (((PM & VP) + VP) ^ VP) | PM | VN;
		HP = VN | ~(D0 | VP);
		HN = D0 & VP;
		if(HP & top) ++D;
		else if(HN & top) --D;
		VP = (HN << 1L) | ~(D0 | ((HP << 1L) | 1L));
		VN = D0 & ((HP << 1L) | 1L);
		// cout << "PM=" << PM << endl;
		// cout << "D0=" << D0 << endl;
		// cout << "HP=" << HP << endl;
		// cout << "HN=" << HN << endl;
		// cout << "VP=" << VP << endl;
		// cout << "VN=" << VN << endl;
		// cout << "D=" << D << endl;
	}
	return D;
}


template<size_t N, typename T, typename TVALUE, typename V>
unsigned int edit_distance_bpv(T &cmap, V const &vec, unsigned int const &tmax, unsigned int const &tlen) {
	int D = tmax * 64 + tlen;
	TVALUE D0, HP, HN, VP = {0}, VN = {0};
	uint64_t top = (1L << (tlen - 1));  // 末尾のvectorに適用
	uint64_t lmb = (1L << 63);
	for(size_t i = 0; i < tmax; ++i) VP[i] = ~0;
	for(size_t i = 0; i < tlen; ++i) VP[tmax] |= (1L << i);
	// cout << "VP0=" << to_str(VP[0]) << endl;
	// cout << "VN0=" << to_str(VN[0]) << endl;
	for(size_t i = 0; i < vec.size(); ++i) {
		TVALUE &PM = cmap[vec[i]];
		for(int r = 0; r <= tmax; ++r) {
			uint64_t X = PM[r];
			if(r > 0 && (HN[r - 1] & lmb)) X |= 1L;
			D0[r] = (((X & VP[r]) + VP[r]) ^ VP[r]) | X | VN[r];
			HP[r] = VN[r] | ~(D0[r] | VP[r]);
			HN[r] = D0[r] & VP[r];
			X = (HP[r] << 1L);
			if(r == 0 || HP[r - 1] & lmb) X |= 1L;
			VP[r] = (HN[r] << 1L) | ~(D0[r] | X);
			if(r > 0 && (HN[r - 1] & lmb)) VP[r] |= 1L;
			VN[r] = D0[r] & X;
			// cout << "r=" << r << endl;
			// cout << "PM(" << vec[i] << ")=" << to_str(PM[r]) << endl;
			// cout << "D0=" << to_str(D0[r]) << endl;
			// cout << "HP=" << to_str(HP[r]) << endl;
			// cout << "HN=" << to_str(HN[r]) << endl;
			// cout << "VP=" << to_str(VP[r]) << endl;
			// cout << "VN=" << to_str(VN[r]) << endl;
		}
		if(HP[tmax] & top) ++D;
		else if(HN[tmax] & top) --D;
		// cout << "D=" << D << endl;
	}
	return D;
}


/// c.f. http://handasse.blogspot.com/2009/04/c_29.html
template<typename T>
unsigned int edit_distance_dp(T const &str1, T const &str2) {
	// vectorより固定長配列の方が速いが、文字列が長い時の保険でのみ呼ばれるのでサイズを決め打ちできない
	vector< vector<uint32_t> > d(str1.size() + 1, vector<uint32_t>(str2.size() + 1));
	for (int i = 0; i < str1.size() + 1; i++) d[i][0] = i;
	for (int i = 0; i < str2.size() + 1; i++) d[0][i] = i;
	for (int i = 1; i < str1.size() + 1; i++) {
		for (int j = 1; j < str2.size() + 1; j++) {
			d[i][j] = min(min(d[i-1][j], d[i][j-1]) + 1, d[i-1][j-1] + (str1[i-1] == str2[j-1] ? 0 : 1));
		}
	}
	return d[str1.size()][str2.size()];
}


template<size_t N, typename T>
unsigned int edit_distance_map_(T const &a, T const &b) {
	typedef map<typename T::value_type, array<uint64_t, N> > cmap_v;
	cmap_v cmap;
	unsigned int tmax = (a.size() - 1) >> 6;
	unsigned int tlen = a.size() - tmax * 64;
	for(size_t i = 0; i < tmax; ++i) {
		for(size_t j = 0; j < 64; ++j) cmap[a[i * 64 + j]][i] |= (1L << j);
	}
	for(size_t i = 0; i < tlen; ++i) cmap[a[tmax * 64 + i]][tmax] |= (1L << i);
	return edit_distance_bpv<N, cmap_v, typename cmap_v::mapped_type, string>(cmap, b, tmax, tlen);
}


/// サイズ１の場合の特殊化
template<typename T>
unsigned int edit_distance_map_1_(T const &a, T const &b) {
	map<char, uint64_t> cmap;
	for(size_t i = 0; i < a.size(); ++i) cmap[a[i]] |= (1L << i);
	return edit_distance_bp<map<char, uint64_t>, string>(cmap, b, a.size());
}


template<size_t N, typename T, size_t M>
unsigned int edit_distance_fixed_(T const &a, T const &b) {
	uint64_t cmap[M][N] = {0};
	unsigned int tmax = (a.size() - 1) >> 6;
	unsigned int tlen = a.size() - tmax * 64;
	for(size_t i = 0; i < tmax; ++i) {
		for(size_t j = 0; j < 64; ++j) cmap[a[i * 64 + j]][i] |= (1L << j);
	}
	for(size_t i = 0; i < tlen; ++i) cmap[a[tmax * 64 + i]][tmax] |= (1L << i);
	return edit_distance_bpv<N, uint64_t[M][N], uint64_t[N], string>(cmap, b, tmax, tlen);
}


/// サイズ１の場合の特殊化
template<typename T, size_t M>
unsigned int edit_distance_fixed_1_(T const &a, T const &b) {
	uint64_t cmap[M] = {0};
	for(size_t i = 0; i < a.size(); ++i) cmap[a[i]] |= (1L << i);
	return edit_distance_bp<uint64_t[M], string>(cmap, b, a.size());
}

template<typename T>
unsigned int edit_distance(T const &a, T const &b, MAPOPTION const &opt) {
	if(a.size() == 0) return b.size();
	else if(b.size() == 0) return a.size();
	// 要素数の大きいほうがa
	T const *ap, *bp;
	if(a.size() < b.size()) ap = &b, bp = &a;
	else ap = &a, bp = &b;
	// 必要な配列サイズを調べる
	size_t vsize = ((ap->size() - 1) >> 6) + 1;  // 64までは1, 128までは2, ...
	// bit-parallelでできそうな限界を超えたら要素数の小さい方をaとする。
	if(vsize > 10) {
		T const *_ = ap;
		ap = bp, bp = _;
		vsize = ((ap->size() - 1) >> 6) + 1;
	}

	switch (vsize) {
		case 1:
			switch(opt) {
				case MAPOPTION::MAP: return edit_distance_map_1_<T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_1_<T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_1_<T, 65536>(*ap, *bp);
			}
		case 2:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<2, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<2, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<2, T, 65536>(*ap, *bp);
			}
		case 3:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<3, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<3, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<3, T, 65536>(*ap, *bp);
			}
		case 4:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<4, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<4, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<4, T, 65536>(*ap, *bp);
			}
		case 5:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<5, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<5, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<5, T, 65536>(*ap, *bp);
			}
		case 6:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<6, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<6, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<6, T, 65536>(*ap, *bp);
			}
		case 7:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<7, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<7, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<7, T, 65536>(*ap, *bp);
			}
		case 8:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<8, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<8, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<8, T, 65536>(*ap, *bp);
			}
		case 9:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<9, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<9, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<9, T, 65536>(*ap, *bp);
			}
		case 10:
			switch (opt) {
				case MAPOPTION::MAP: return edit_distance_map_<10, T>(*ap, *bp);
				case MAPOPTION::N1: return edit_distance_fixed_<10, T, 256>(*ap, *bp);
				case MAPOPTION::N2: return edit_distance_fixed_<10, T, 65536>(*ap, *bp);
			}
		default:
			return edit_distance_dp<T>(*ap, *bp);  // dynamic programmingに任せる
	}
}

unsigned int edit_distance_c(const char *a, const size_t asize, const char *b, const size_t bsize) {
	const string sa(a, asize);
	const string sb(b, bsize);
	return edit_distance<string>(sa, sb, MAPOPTION::MAP);
}

void edit_distance_matrix(char **words, const size_t *sizes, size_t count) {
	unsigned int dist = 0;
	vector<string> cpp_words(count);
	cpp_words.resize(count);
	for (size_t i = 0; i < count; ++i) {
		cpp_words[i] = string(words[i], sizes[i]);
	}

	for (size_t i = 0; i < count; ++i) {
		const string word = cpp_words[i];
		cout << "." << flush;
		for (size_t j = 0; j < count; ++j) {
			dist += edit_distance<string>(word, cpp_words[j], MAPOPTION::MAP);
		}
	}
	cout << dist << endl;
}
