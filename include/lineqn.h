#ifndef LINEQN_H
#define LINEQN_H
/*
Szymon Rusinkiewicz
Princeton University

lineqn.h
Solution of systems of linear equations and eigenvalue decomposition.
*/


#include "mathutil.h"


// GCC generates (painstakingly-verified-to-be) bogus warnings in this file
#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 406)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
# pragma GCC diagnostic ignored "-Warray-bounds"
#endif


namespace trimesh {

// LU decomposition with partial pivoting.
// A overwritten with L and U in the usual way, while ind gets the permutation.
// Returns true iff succeeded (i.e., A is nonsingular).
// Patterned after Numerical Recipes.
template <class T, int N>
static inline bool ludcmp(T (&A)[N][N], int ind[N])
{
	using namespace ::std;
	T pivot[N];

	for (int i = 0; i < N; i++) {
		T big = 0;
		for (int j = 0; j < N; j++) {
			T tmp = abs(A[i][j]);
			if (tmp > big)
				big = tmp;
		}
		if (big == 0)
			return false;
		pivot[i] = 1 / big;
	}
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < j; i++) {
			T sum = A[i][j];
			for (int k = 0; k < i; k++)
				sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
		}
		T big = 0;
		int imax = j;
		for (int i = j; i < N; i++) {
			T sum = A[i][j];
			for (int k = 0; k < j; k++)
				sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
			T tmp = pivot[i] * abs(sum);
			if (tmp > big) {
				big = tmp;
				imax = i;
			}
		}
		if (imax != j) {
			for (int k = 0; k < N; k++)
				swap(A[imax][k], A[j][k]);
			pivot[imax] = pivot[j];
		}
		ind[j] = imax;
		if (unlikely(A[j][j] == 0))
			return false;
		if (j != N - 1) {
			T tmp = 1 / A[j][j];
			for (int i = j + 1; i < N; i++)
				A[i][j] *= tmp;
		}
	}
	return true;
}


// Solve Ax=b after ludcmp via backsubstitution.  x overwrites b.
template <class T, int N>
static inline void lubksb(const T (&A)[N][N], const int ind[N], T b[N])
{
	for (int i = 0; i < N; i++) {
		int ip = ind[i];
		T sum = b[ip];
		b[ip] = b[i];
		for (int j = 0; j < i; j++)
			sum -= A[i][j] * b[j];
		b[i] = sum;
	}
	for (int i = N - 1; i >= 0; i--) {
		T sum = b[i];
		for (int j = i + 1; j < N; j++)
			sum -= A[i][j] * b[j];
		b[i] = sum / A[i][i];
	}
}


// Solve Ax=b after ludcmp via backsubstitution, with x separate from b.
template <class T, int N>
static inline void lubksb(const T (&A)[N][N],
                          const int ind[N],
                          const T b[N],
                          T x[N])
{
	if (&x[0] != &b[0]) {
		for (int i = 0; i < N; i++)
			x[i] = b[i];
	}
	lubksb(A, ind, x);
}


// LDL^T decomposition of a symmetric positive definite matrix (and some
// other symmetric matrices, but fragile since we don't do pivoting).
// Like Cholesky, but no square roots, which is important for small N.
// Reads diagonal and upper triangle of matrix A.
// On output, lower triangle of A holds LD, while rdiag holds D^-1.
// Algorithm from Golub and van Loan.
template <class T, int N>
static inline bool ldltdc(T (&A)[N][N], T rdiag[N])
{
	// Special case for small N
	if (N < 1) {
		return false;
	} else if (N <= 3) {
		T d0 = A[0][0];
		rdiag[0] = 1 / d0;
		if (N == 1)
			return (d0 != 0);
		A[1][0] = A[0][1];
		T l10 = rdiag[0] * A[1][0];
		T d1 = A[1][1] - l10 * A[1][0];
		rdiag[1] = 1 / d1;
		if (N == 2)
			return (d0 != 0 && d1 != 0);
		T d2 = A[2][2] - rdiag[0] * sqr(A[2][0]) - rdiag[1] * sqr(A[2][1]);
		rdiag[2] = 1 / d2;
		A[2][0] = A[0][2];
		A[2][1] = A[1][2] - l10 * A[2][0];
		return (d0 != 0 && d1 != 0 && d2 != 0);
	}

	T v[N-1];
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < i; k++)
			v[k] = A[i][k] * rdiag[k];
		for (int j = i; j < N; j++) {
			T sum = A[i][j];
			for (int k = 0; k < i; k++)
				sum -= v[k] * A[j][k];
			if (i == j) {
				if (unlikely(sum == 0))
					return false;
				rdiag[i] = 1 / sum;
			} else {
				A[j][i] = sum;
			}
		}
	}

	return true;
}


// Solve Ax=b after ldltdc.  x is allowed to be the same as b.
template <class T, int N>
static inline void ldltsl(const T (&A)[N][N],
                          const T rdiag[N],
                          const T b[N],
                          T x[N])
{
	for (int i = 0; i < N; i++) {
		T sum = b[i];
		for (int k = 0; k < i; k++)
			sum -= A[i][k] * x[k];
		x[i] = sum * rdiag[i];
	}
	for (int i = N - 1; i >= 0; i--) {
		T sum = 0;
		for (int k = i + 1; k < N; k++)
			sum += A[k][i] * x[k];
		x[i] -= sum * rdiag[i];
	}
}


// Convenience function when output x overwrites input b.
template <class T, int N>
static inline void ldltsl(const T (&A)[N][N], const T rdiag[N], T b[N])
{
	ldltsl(A, rdiag, b, b);
}


// Eigenvector decomposition for real, symmetric matrices.
// Entries of d are eigenvalues, sorted smallest to largest.
// A changed in-place to have its columns hold the corresponding eigenvectors.
// Note that A must be completely filled in on input.
// Bowdler et al. algorithm, via EISPACK / JAMA / TNT.
template <class T, int N>
static inline void eigdc(T (&A)[N][N], T d[N])
{
	using namespace ::std;
	T e[N];

	// Householder
	for (int j = 0; j < N; j++) {
		d[j] = A[N-1][j];
		e[j] = 0;
	}
	for (int i = N - 1; i > 0; i--) {
		T scale = 0;
		for (int k = 0; k < i; k++)
			scale += abs(d[k]);
		if (scale == 0) {
			e[i] = d[i-1];
			for (int j = 0; j < i; j++) {
				d[j] = A[i-1][j];
				A[i][j] = A[j][i] = 0;
			}
			d[i] = 0;
		} else {
			T h(0);
			T invscale = 1 / scale;
			for (int k = 0; k < i; k++) {
				d[k] *= invscale;
				h += sqr(d[k]);
			}
			T f = d[i-1];
			T g = (f > 0) ? -sqrt(h) : sqrt(h);
			e[i] = scale * g;
			h -= f * g;
			d[i-1] = f - g;
			for (int j = 0; j < i; j++)
				e[j] = 0;
			for (int j = 0; j < i; j++) {
				f = d[j];
				A[j][i] = f;
				g = e[j] + f * A[j][j];
				for (int k = j + 1; k < i; k++) {
					g += A[k][j] * d[k];
					e[k] += A[k][j] * f;
				}
				e[j] = g;
			}
			f = 0;
			T invh = 1 / h;
			for (int j = 0; j < i; j++) {
				e[j] *= invh;
				f += e[j] * d[j];
			}
			T hh = f / (h + h);
			for (int j = 0; j < i; j++)
				e[j] -= hh * d[j];
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k < i; k++)
					A[k][j] -= f * e[k] + g * d[k];
				d[j] = A[i-1][j];
				A[i][j] = 0;
			}
			d[i] = h;
		}
	}

	for (int i = 0; i < N - 1; i++) {
		A[N-1][i] = A[i][i];
		A[i][i] = 1;
		T h = d[i+1];
		if (h != 0) {
			T invh = 1 / h;
			for (int k = 0; k <= i; k++)
				d[k] = A[k][i+1] * invh;
			for (int j = 0; j <= i; j++) {
				T g = 0;
				for (int k = 0; k <= i; k++)
					g += A[k][i+1] * A[k][j];
				for (int k = 0; k <= i; k++)
					A[k][j] -= g * d[k];
			}
		}
		for (int k = 0; k <= i; k++)
			A[k][i+1] = 0;
	}
	for (int j = 0; j < N; j++) {
		d[j] = A[N-1][j];
		A[N-1][j] = 0;
	}
	A[N-1][N-1] = 1;

	// QL
	for (int i = 1; i < N; i++)
		e[i-1] = e[i];
	e[N-1] = 0;
	T f = 0, tmp = 0;
	const T eps = numeric_limits<T>::epsilon();
	for (int l = 0; l < N; l++) {
		tmp = max(tmp, abs(d[l]) + abs(e[l]));
		int m = l;
		while (m < N) {
			if (abs(e[m]) <= eps * tmp)
				break;
			m++;
		}
		if (m > l) {
			do {
				T g = d[l];
				T p = (d[l+1] - g) / (e[l] + e[l]);
				T r = hypot(p, T(1));
				if (p < 0)
					r = -r;
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				T dl1 = d[l+1];
				T h = g - d[l];
				for (int i = l + 2; i < N; i++)
					d[i] -= h;
				f += h;
				p = d[m];
				T c = 1, c2 = 1, c3 = 1;
				T el1 = e[l+1], s = 0, s2 = 0;
				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot(p, e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);
					for (int k = 0; k < N; k++) {
						h = A[k][i+1];
						A[k][i+1] = s * A[k][i] + c * h;
						A[k][i] = c * A[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;
			} while (abs(e[l]) > eps * tmp);
		}
		d[l] += f;
		e[l] = 0;
	}

	// Sort
	for (int i = 0; i < N - 1; i++) {
		int k = i;
		T p = d[i];
		for (int j = i+1; j < N; j++) {
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}
		if (k == i)
			continue;
		d[k] = d[i];
		d[i] = p;
		for (int j = 0; j < N; j++)
			swap(A[j][i], A[j][k]);
	}
}


// x <- A * d * A' * b.  x is allowed to be the same as b.
template <class T, int N>
static inline void eigmult(const T (&A)[N][N],
                           const T d[N],
                           const T b[N],
                           T x[N])
{
	T e[N];
	for (int i = 0; i < N; i++) {
		T sum = 0;
		for (int j = 0; j < N; j++)
			sum += A[j][i] * b[j];
		e[i] = sum * d[i];
	}
	for (int i = 0; i < N; i++) {
		T sum = 0;
		for (int j = 0; j < N; j++)
			sum += A[i][j] * e[j];
		x[i] = sum;
	}
}


// Convenience function when output x overwrites input b.
template <class T, int N>
static inline void eigmult(const T (&A)[N][N], const T d[N], T b[N])
{
	eigmult<T,N>(A, d, b, b);
}


// Singular value decomposition A -> U * diag(s) * V'
// On output, U overwrites A and entries of s are sorted largest to smallest.
// Golub-Reinsch algorithm, via LINPACK / JAMA / TNT.
template <class T, int M, int N>
static inline void svd(T (&A)[M][N], T s[N], T (&V)[N][N])
{
	using namespace ::std;
	T e[N], work[M];

	const int ncu = min(M, N);
	const int nct = min(M - 1, N);
	const int nrt = max(min(M, N - 2), 0);
	const int lu = max(nct, nrt);

	// Reduce A to bidiagonal form
	for (int k = 0; k < lu; k++) {
		if (k < nct) {
			s[k] = 0;
			for (int i = k; i < M; i++)
				s[k] = hypot(s[k], A[i][k]);
			if (s[k] != 0) {
				if (A[k][k] < 0)
					s[k] = -s[k];
				for (int i = k; i < M; i++)
					A[i][k] /= s[k];
				A[k][k] += 1;
			}
			s[k] = -s[k];
		}
		for (int j = k + 1; j < N; j++) {
			if (k < nct && s[k] != 0) {
				T t = 0;
				for (int i = k; i < M; i++)
					t += A[i][k] * A[i][j];
				t = -t / A[k][k];
				for (int i = k; i < M; i++)
					A[i][j] += t * A[i][k];
			}
			e[j] = A[k][j];
		}
		if (k < nrt) {
			e[k] = 0;
			for (int i = k + 1; i < N; i++)
				e[k] = hypot(e[k], e[i]);
			if (e[k] != 0) {
				if (e[k+1] < 0)
					e[k] = -e[k];
				for (int i = k + 1; i < N; i++)
					e[i] /= e[k];
				e[k+1] += 1;
			}
			e[k] = -e[k];
			if (k+1 < M && e[k] != 0) {
				for (int i = k + 1; i < M; i++)
					work[i] = 0;
				for (int j = k + 1; j < N; j++) {
					for (int i = k + 1; i < M; i++)
						work[i] += e[j] * A[i][j];
				}
				for (int j = k + 1; j < N; j++) {
					T t = -e[j] / e[k+1];
					for (int i = k + 1; i < M; i++)
						A[i][j] += t * work[i];
				}
			}
			for (int i = k + 1; i < N; i++)
				V[i][k] = e[i];
		}
	}

	int p = min(N, M + 1);
	if (nct < N)
		s[nct] = A[nct][nct];
	if (p > M)
		s[p-1] = 0;
	if (nrt + 1 < p)
		e[nrt] = A[nrt][p-1];
	e[p-1] = 0;

	// Replace A with U
	for (int i = 0; i < ncu; i++)
		for (int j = i + 1; j < N; j++)
			A[i][j] = 0;
	if (nct < ncu)
		A[nct][nct] = 1;
	for (int k = nct - 1; k >= 0; k--) {
		if (s[k] != 0) {
			for (int j = k + 1; j < ncu; j++) {
				T t = 0;
				for (int i = k; i < M; i++)
					t += A[i][k] * A[i][j];
				t = -t / A[k][k];
				for (int i = k; i < M; i++)
					A[i][j] += t * A[i][k];
			}
			for (int i = k; i < M; i++)
				A[i][k] = -A[i][k];
			A[k][k] += 1;
			for (int i = 0; i < k - 1; i++)
				A[i][k] = 0;
		} else {
			for (int i = 0; i < M; i++)
				A[i][k] = 0;
			A[k][k] = 1;
		}
	}

	// Generate V
	for (int k = N - 1; k >= 0; k--) {
		if (k < nrt && e[k] != 0) {
			for (int j = k + 1; j < N; j++) {
				T t = 0;
				for (int i = k + 1; i < N; i++)
					t += V[i][k] * V[i][j];
				t = -t / V[k+1][k];
				for (int i = k + 1; i < N; i++)
					V[i][j] += t * V[i][k];
			}
		}
		for (int i = 0; i < N; i++)
			V[i][k] = 0;
		V[k][k] = 1;
	}


	// Main iteration loop for the singular values.
	const T eps = numeric_limits<T>::epsilon();
	const T tiny = numeric_limits<T>::min() / eps;
	int pp = p - 1;
	while (p > 0) {
		int k = p - 2;
		for ( ; k >= 0; k--) {
			if (abs(e[k]) <=
			    tiny + eps * (abs(s[k]) + abs(s[k+1]))) {
				e[k] = 0;
				break;
			}
		}

		int kase = 4;
		if (k != p - 2) {
			int ks = p - 1;
			for ( ; ks > k; ks--) {
				T t = (ks != p ? abs(e[ks]) : T(0)) +
				    (ks != k+1 ? abs(e[ks-1]) : T(0));
				if (abs(s[ks]) <= tiny + eps * t) {
					s[ks] = 0;
					break;
				}
			}
			if (ks == k) {
				kase = 3;
			} else if (ks == p-1) {
				kase = 1;
			} else {
				kase = 2;
				k = ks;
			}
		}
		k++;

		switch (kase) {
		case 1: {
			T f = e[p-2];
			e[p-2] = 0;
			for (int j = p - 2; j >= k; j--) {
				T t = hypot(s[j], f);
				T cs = s[j] / t;
				T sn = f / t;
				s[j] = t;
				if (j != k) {
					f = -sn * e[j-1];
					e[j-1] = cs * e[j-1];
				}
				for (int i = 0; i < N; i++) {
					t = cs * V[i][j] + sn * V[i][p-1];
					V[i][p-1] =
					   -sn * V[i][j] + cs * V[i][p-1];
					V[i][j] = t;
				}
			}
		}
		break;

		case 2: {
			T f = e[k-1];
			e[k-1] = 0;
			for (int j = k; j < p; j++) {
				T t = hypot(s[j], f);
				T cs = s[j] / t;
				T sn = f / t;
				s[j] = t;
				f = -sn * e[j];
				e[j] = cs * e[j];
				for (int i = 0; i < M; i++) {
					t = cs * A[i][j] + sn * A[i][k-1];
					A[i][k-1] =
					   -sn * A[i][j] + cs * A[i][k-1];
					A[i][j] = t;
				}
			}
		}
		break;

		case 3: {
			T scale = 1 / max(max(max(max(
				abs(s[p-1]), abs(s[p-2])), abs(e[p-2])),
				abs(s[k])), abs(e[k]));
			T sp = s[p-1] * scale;
			T spm1 = s[p-2] * scale;
			T epm1 = e[p-2] * scale;
			T sk = s[k] * scale;
			T ek = e[k] * scale;
			T b = 0.5 * ((spm1 + sp) * (spm1 - sp) + epm1 * epm1);
			T c = sqr(sp * epm1);
			T shift = 0;
			if ((b != 0) || (c != 0)) {
				shift = sqrt(b * b + c);
				if (b < 0)
					shift = -shift;
				shift = c / (b + shift);
			}
			T f = (sk + sp) * (sk - sp) + shift;
			T g = sk * ek;

			for (int j = k; j < p - 1; j++) {
				T t = hypot(f, g);
				T cs = f / t;
				T sn = g / t;
				if (j != k)
					e[j-1] = t;
				f = cs * s[j] + sn * e[j];
				e[j] = cs * e[j] - sn * s[j];
				g = sn * s[j+1];
				s[j+1] = cs * s[j+1];
				for (int i = 0; i < N; i++) {
					t = cs * V[i][j] + sn * V[i][j+1];
					V[i][j+1] =
					   -sn * V[i][j] + cs * V[i][j+1];
					V[i][j] = t;
				}
				t = hypot(f, g);
				cs = f / t;
				sn = g / t;
				s[j] = t;
				f = cs * e[j] + sn * s[j+1];
				s[j+1] =
				   -sn * e[j] + cs * s[j+1];
				g = sn * e[j+1];
				e[j+1] = cs * e[j+1];
				if (j >= M - 1)
					continue;
				for (int i = 0; i < M; i++) {
					t = cs * A[i][j] + sn * A[i][j+1];
					A[i][j+1] =
					   -sn * A[i][j] + cs * A[i][j+1];
					A[i][j] = t;
				}
			}
			e[p-2] = f;
		}
		break;

		// Convergence.
		case 4: {
			// Make the singular values positive.
			if (s[k] <= 0) {
				s[k] = -s[k];
				for (int i = 0; i < N; i++)
					V[i][k] = -V[i][k];
			}

			// Order the singular values.
			while (k < pp) {
				if (s[k] >= s[k+1])
					break;
				T t = s[k];
				s[k] = s[k+1];
				s[k+1] = t;
				if (k < N - 1) {
					for (int i = 0; i < N; i++) {
						t = V[i][k+1];
						V[i][k+1] = V[i][k];
						V[i][k] = t;
					}
				}
				if (k < M - 1) {
					for (int i = 0; i < M; i++) {
						t = A[i][k+1];
						A[i][k+1] = A[i][k];
						A[i][k] = t;
					}
				}
				k++;
			}
			p--;
		}
		break;

		} // switch
	} // while (p > 0)
}

} // namespace trimesh

#if defined(__GNUC__) && (__GNUC__ * 100 + __GNUC_MINOR__ >= 406)
# pragma GCC diagnostic pop
#endif

#endif
