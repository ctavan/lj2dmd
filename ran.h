#ifndef RAN_H
#define RAN_H
/**
 * Press, Teukolsky, Vetterling, Flannery
 * NUMERICAL RECIPES - The Art of Scientific Computing - Third Edition (2007)
 * p. 342 / 343
 */
typedef int Int; // 32 bit integer
typedef unsigned int Uint;

typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;

typedef double Doub; // default floating type
typedef long double Ldoub;

struct Ran {
	Ullong u,v,w;
	Ran(Ullong j) : v(4101842887655102017LL), w(1) {
		u = j ^ v; int64();
		v = u; int64();
		w = v; int64();
	}
	inline Ullong int64() {
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
		w = 4294957665U*(w & 0xffffffff) + (w >> 32);
		Ullong x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
		return (x + v) ^ w;
	}
	inline Doub doub() { return 5.42101086242752217E-20 * int64(); }
	inline Uint int32() { return (Uint)int64(); }
};

#endif
