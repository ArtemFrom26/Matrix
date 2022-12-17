#include "biginteger.h"
#include <array>
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t N, size_t K, bool flag>
struct IsPrimeHelper {
	static const bool value = (N % K != 0) && IsPrimeHelper<N, K + 1, K * K <= N>::value;
};
template <size_t N, size_t K>
struct IsPrimeHelper<N, K, false> {
	static const bool value = true;
};
template <size_t N>
struct IsPrimeHelper<N, 1, true> {
	static const bool value = IsPrimeHelper<N, 2, 4 <= N>::value;
};
template <size_t N>
struct IsPrime {
	static const bool value = IsPrimeHelper<N, 1, true>::value;
};

template <size_t N, size_t L>
struct MaxHelper {
	static const size_t value = (N > L ? N : L);
};
template <size_t M, size_t N, size_t L>
struct Max {
	static const size_t value = (M > MaxHelper<N, L>::value ? M : MaxHelper<N, L>::value);
};

template <size_t N, size_t current, bool flag>
struct Helper;
template <size_t N, size_t current>
struct Helper<N, current, false> {
	static const size_t value = Helper<N, current * 2, N <= 2 * current>::value;
};
template <size_t N, size_t current>
struct Helper<N, current, true> {
	static const size_t value = current;
};
template <size_t N>
struct ClosestTwoDegree{
	static const size_t value = Helper<N, 1, false>::value;
};
template <>
struct ClosestTwoDegree<1> {
	static const size_t value = 1;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t N>
class Residue {
  private:
	long long value;

  public:
  	Residue();
  	Residue(long long n);
  	explicit operator int() const;
  	
  	long long getValue() const;
  	long long& getValue();
  	std::string toString() const;
  	
	Residue& operator+=(const Residue& right);
    	Residue& operator-=(const Residue& right);
	Residue& operator*=(const Residue& right);
	Residue& operator/=(const Residue& right);
	Residue operator-() const;
};

template <size_t N>
Residue<N> operator+(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
Residue<N> operator-(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
Residue<N> operator*(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
Residue<N> operator/(const Residue<N>& left, const Residue<N>& right);

template <size_t N>
bool operator==(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
bool operator!=(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
bool operator<(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
bool operator<=(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
bool operator>(const Residue<N>& left, const Residue<N>& right);
template <size_t N>
bool operator>=(const Residue<N>& left, const Residue<N>& right);

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& num);
template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& num);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t N>
Residue<N>::Residue() : value(0) {
}	
template <size_t N>
Residue<N>::Residue(long long n) {
	if (n < 0) {
		n += ((-n)/N + 1) * N;
	}

	value = n % N;
}
template <size_t N>
Residue<N>::operator int() const {
	return static_cast<int>(value);
}
	
template <size_t N>
long long Residue<N>::getValue() const {
	return value;
}
template <size_t N>
long long& Residue<N>::getValue() {
	return value;
}
template <size_t N>
std::string Residue<N>::toString() const {
	return std::to_string(value);
}
	
template <size_t N>
Residue<N>& Residue<N>::operator+=(const Residue<N>& right) {
	value = (value + right.value) % N;
	return *this;
}
template <size_t N>
Residue<N>& Residue<N>::operator-=(const Residue<N>& right) {
	value -= right.value;
	
	if (value < 0) {
		value += N;
	}

	return *this;
}
template <size_t N>
Residue<N>& Residue<N>::operator*=(const Residue<N>& right) {
	value *= right.value;
	value %= N;

	return *this;
}
template <size_t N>
Residue<N>& Residue<N>::operator/=(const Residue<N>& right) {
	static_assert(IsPrime<N>::value);
		
	long long val = right.getValue();
	long long inverse = 1;
	size_t n = N - 2;

	while (n != 0) {
		if (n & 1) {
			inverse *= val;
			inverse %= N;
		}

		val *= val;
		val %= N;
		n >>= 1;
	}
		
	value *= inverse;
	value %= N;

	return *this;
}
template <size_t N>
Residue<N> Residue<N>::operator-() const {
	return {(long long)(N - value)};
}

template <size_t N>
Residue<N> operator+(const Residue<N>& left, const Residue<N>& right) {
 	Residue<N> tmp = left;
 	tmp += right;
	return tmp;
}
template <size_t N>
Residue<N> operator-(const Residue<N>& left, const Residue<N>& right) {
  	Residue<N> tmp = left;
  	tmp -= right;
  	return tmp;
}
template <size_t N>
Residue<N> operator*(const Residue<N>& left, const Residue<N>& right) {
  	Residue<N> tmp = left;
  	tmp *= right;
  	return tmp;
}
template <size_t N>
Residue<N> operator/(const Residue<N>& left, const Residue<N>& right) {
  	Residue<N> tmp = left;
  	tmp /= right;
  	return tmp;
}

template <size_t N>
bool operator==(const Residue<N>& left, const Residue<N>& right) {
  	return left.getValue() == right.getValue();
}
template <size_t N>
bool operator!=(const Residue<N>& left, const Residue<N>& right) {
  	return left.getValue() != right.getValue();
}
template <size_t N>
bool operator<(const Residue<N>& left, const Residue<N>& right) {
  	return left.getValue() < right.getValue();
}
template <size_t N>
bool operator<=(const Residue<N>& left, const Residue<N>& right) {
  	return left.getValue() <= right.getValue();
}
template <size_t N>
bool operator>(const Residue<N>& left, const Residue<N>& right) {
  	return left.getValue() > right.getValue();
}
template <size_t N>
bool operator>=(const Residue<N>& left, const Residue<N>& right) {
 	return left.getValue() >= right.getValue();
}

template <size_t N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& num) {
	out << num.toString();
	return out;
}
template <size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& num) {
	long long n;
	in >> n;
		
	if (n < 0) {
		n += ((-n)/N + 1) * N;
	}
	n %= N;
	num.getValue() = n;
	
	return in;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t M, size_t N, typename Field = Rational>
class Matrix {
  private:
	std::array<std::array<Field, N>, M> body; //массив строк	

  public:
	Matrix();
	Matrix(std::initializer_list<std::array<Field, N>> list);
	
	const std::array<Field, N>& operator[](size_t idx) const;
	std::array<Field, N>& operator[](size_t idx);
	std::vector<Field> getRow(size_t idx);
	std::vector<Field> getColumn(size_t idx);
	
	Matrix<M, N, Field>& operator*=(const Field& coefficient);
	Matrix<M, N, Field>& operator+=(const Matrix<M, N, Field>& right);
	Matrix<M, N, Field>& operator-=(const Matrix<M, N, Field>& right);
	template <size_t L>
	Matrix<M, L, Field>& operator*=(const Matrix<N, L, Field>& right);
	
	size_t rank() const;
	Field trace() const;
	Matrix<N, M, Field> transposed() const;
	Field det() const;
	Matrix<M, N, Field> inverted() const;
	void invert();
};
template<size_t N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;

template<size_t N, typename Field>
SquareMatrix<N, Field> multiply(const SquareMatrix<N, Field>& left, const SquareMatrix<N, Field>& right);

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right);
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right);
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& matrix, const Field& coefficient);
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& coefficient, const Matrix<M, N, Field>& matrix);
template <size_t M, size_t N, size_t L, typename Field>
Matrix<M, L, Field> operator*(const Matrix<M, N, Field>& left, const Matrix<N, L, Field>& right);

template <size_t M, size_t N, typename Field, size_t K, size_t L>
bool operator==(const Matrix<M, N, Field>& left, const Matrix<K, L, Field>& right);
template <size_t M, size_t N, typename Field, size_t K, size_t L>
bool operator!=(const Matrix<M, N, Field>& left, const Matrix<K, L, Field>& right);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <size_t M, size_t N, typename Field>	
Matrix<M, N, Field>::Matrix() {
	for (size_t i = 0; i < M; ++i) {
		body[i].fill(0);
	}
}
template <size_t M, size_t N, typename Field>	
Matrix<M, N, Field>::Matrix(std::initializer_list<std::array<Field, N>> list) {
	auto it = list.begin();
	for (size_t i = 0; i < M; ++i, ++it) {
		for (size_t j = 0; j < N; ++j) {
			body[i][j] = (*it)[j];
		}
	}
}
    	
template <size_t M, size_t N, typename Field>
const std::array<Field, N>& Matrix<M, N, Field>::operator[](size_t idx) const {
	return body[idx];
}
template <size_t M, size_t N, typename Field>
std::array<Field, N>& Matrix<M, N, Field>::operator[](size_t idx) {
	return body[idx];
}
template <size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getRow(size_t idx) {
	std::vector<Field> result(N);
	for (size_t i = 0; i < N; ++i) {
		result[i] = body[idx][i];
	}

	return result;
} 
template <size_t M, size_t N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(size_t idx) {
	std::vector<Field> result(M);
	for (size_t i = 0; i < M; ++i) {
		result[i] = body[i][idx];
	}

	return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator*=(const Field& coefficient) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			body[i][j] *= coefficient;
		}
	}
	
	return *this;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field>& right) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			body[i][j] += right[i][j];
		}
	}

	return *this;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field>& Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field>& right) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			body[i][j] -= right[i][j];
		}
	}
		
	return *this;
}
template <size_t M, size_t N, typename Field>
template <size_t L>
Matrix<M, L, Field>& Matrix<M, N, Field>::operator*=(const Matrix<N, L, Field>& right) {
	if (M == N && L == N) {
		*this = multiply(*this, right);
		return *this;
	}
	
	*this = *this * right;
	return *this;
}
template <size_t M, size_t N, typename Field>
size_t Matrix<M, N, Field>::rank() const {
	size_t rank = 0;
	std::vector<bool> line_used(M, false);
	std::array<std::array<Field, N>, M> temp = body;

	for (size_t i = 0; i < N; ++i) {
		size_t j;
		for (j = 0; j < M; ++j) {
			if (!line_used[j] && temp[j][i] != Field(0)) {
				break;
			}
		}
		
		if (j != M) {
			++rank;
			line_used[j] = true;
			for (size_t p = i + 1; p < N; ++p) {
				temp[j][p] /= temp[j][i];
			}
			for (size_t k = 0; k < M; ++k) {
				if (k != j && temp[j][i] != Field(0)) {
					for (size_t p = i + 1; p < N; ++p) {
						temp[k][p] -= temp[j][p] * temp[k][i];
					}
				}
			}
		}
	}

	return rank;
}
template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::trace() const {
	static_assert(M == N);
	
	Field result = 0;
	for (size_t i = 0; i < N; ++i) {
		result += body[i][i];	
	}
	
	return result;
}
template <size_t M, size_t N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
	Matrix<N, M, Field> result;

	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			result[j][i] = body[i][j];
		}
	}

	return result;
}
template <size_t M, size_t N, typename Field>
Field Matrix<M, N, Field>::det() const {
	static_assert(M == N);
	
	Field result = 1;
	std::array<std::array<Field, N>, M> temp = body;
	for (size_t i = 0; i < N; ++i) {
		size_t k = i;
		for (size_t j = i + 1; j < N; ++j) {
			if (temp[j][i] > temp[k][i]) {
				k = j;
			}
		}
		if (temp[k][i] == Field(0)) {
			result = 0;
			break;
		}
		std::swap(temp[i], temp[k]);
		if (i != k) {
			result = -result;
		}
		result *= temp[i][i];
		for (size_t j = i + 1; j < N; ++j) {
			temp[i][j] /= temp[i][i];
		}
		for (size_t j = 0; j < N; ++j) {
			if (j != i && temp[j][i] != Field(0)) {
				for (size_t k = i + 1; k < N; ++k) {
					temp[j][k] -= temp[i][k] * temp[j][i];
				}
			}
		}
	}

	return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
	static_assert(M == N);
	
	Matrix<N, N, Field> temp = *this;
	temp.invert();
	
	return temp;
}
template <size_t M, size_t N, typename Field>
void Matrix<M, N, Field>::invert() {
	static_assert(M == N);
	if (det() == Field(0)) {
		std::cout << "determinant is 0!\n";
		return;
	}
	
    SquareMatrix<N, Field> E;
    for (size_t i = 0; i < N; ++i) {
        E[i][i] = 1;
	}
 	
	Field temp;
    for (size_t k = 0; k < N; ++k) {
        temp = body[k][k];
 
        for (size_t j = 0; j < N; ++j) {
            body[k][j] /= temp;
            E[k][j] /= temp;
        }
 
        for (size_t i = k + 1; i < N; ++i) {
            temp = body[i][k];
 
            for (size_t j = 0; j < N; ++j) {
                body[i][j] -= body[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }
 
    for (size_t k = N - 1; k > 0; --k) {
        for (size_t i = k - 1; i >= 0; --i) {
            temp = body[i][k];
 
            for (size_t j = 0; j < N; ++j) {
                body[i][j] -= body[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
			if (i == 0) {
				break;
			}
        }
    }
 
    *this = E;
}

template <size_t N, typename Field>
SquareMatrix<N, Field> multiply(const SquareMatrix<N, Field>& left, const SquareMatrix<N, Field>& right) {
	if (N <= 64) {
		SquareMatrix<N, Field> result;
		
		for (size_t i = 0; i < N; ++i) {
    			for (size_t j = 0; j < N; ++j) {
					for (size_t k = 0; k < N; ++k) {
	    				result[i][j] += left[i][k] * right[k][j];
					}
    			}
		}

		return result;
	}
					     
	SquareMatrix<N / 2, Field> A11, A12, A21, A22, B11, B12, B21, B22;
	for (size_t i = 0; i < N / 2; ++i) {
		for (size_t j = 0; j < N / 2; ++j) {
			A11[i][j] = left[i][j];
			B11[i][j] = right[i][j];
		}
	}
	for (size_t i = 0; i < N / 2; ++i) {
		for (size_t j = N / 2; j < N; ++j) {
			A12[i][j - N / 2] = left[i][j];
			B12[i][j - N / 2] = right[i][j];
		}
	}
	for (size_t i = N / 2; i < N; ++i) {
		for (size_t j = 0; j < N / 2; ++j) {
			A21[i - N / 2][j] = left[i][j];
			B21[i - N / 2][j] = right[i][j];
		}
	}
	for (size_t i = N / 2; i < N; ++i) {
		for (size_t j = N / 2; j < N; ++j) {
			A22[i - N / 2][j - N / 2] = left[i][j];
			B22[i - N / 2][j - N / 2] = right[i][j];
		}
	}

	SquareMatrix<N / 2, Field> P1 = multiply(A11 + A22, B11 + B22),
				   			   P2 = multiply(A21 + A22, B11),
							   P3 = multiply(A11, B12 - B22),
				   			   P4 = multiply(A22, B21 - B11),
				   			   P5 = multiply(A11 + A12, B22),
				   			   P6 = multiply(A21 - A11, B11 + B12),
				   			   P7 = multiply(A12 - A22, B21 + B22);
	SquareMatrix<N / 2, Field> C11 = P1 + P4 - P5 + P7,
				   			   C12 = P3 + P5,
				  			   C21 = P2 + P4,
				   			   C22 = P1 - P2 + P3 + P6;

	SquareMatrix<N, Field> result;
	for (size_t i = 0; i < N / 2; ++i) {
		for (size_t j = 0; j < N / 2; ++j) {
			result[i][j] = C11[i][j];
		}
	}
	for (size_t i = 0; i < N / 2; ++i) {
		for (size_t j = N / 2; j < N; ++j) {
			result[i][j] = C12[i][j - N / 2];
		}
	}
	for (size_t i = N / 2; i < N; ++i) {
		for (size_t j = 0; j < N / 2; ++j) {
			result[i][j] = C21[i - N / 2][j];
		}
	}
	for (size_t i = N / 2; i < N; ++i) {
		for (size_t j = N / 2; j < N; ++j) {
			result[i][j] = C22[i - N / 2][j - N / 2];
		}
	}
	
	return result;
}

template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right) {
	Matrix<M, N, Field> result(left);
	result += right;
	return result;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field>& left, const Matrix<M, N, Field>& right) {
	Matrix<M, N, Field> result(left);
	result -= right;
	return result;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& matrix, const Field& coefficient) {
	Matrix<M, N, Field> result(matrix);
	result *= coefficient;
	return result;
}
template <size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& coefficient, const Matrix<M, N, Field>& matrix) {
	Matrix<M, N, Field> result(matrix);
	result *= coefficient;
	return result;
}
template <size_t M, size_t N,  size_t L, typename Field>
Matrix<M, L, Field> operator*(const Matrix<M, N, Field>& left, const Matrix<N, L, Field>& right) {
	if (Max<M, N, L>::value <= 64) {
		Matrix<M, L, Field> result;
		
		for (size_t i = 0; i < M; ++i) {
    			for (size_t j = 0; j < L; ++j) {
					for (size_t k = 0; k < N; ++k) {
	    				result[i][j] += left[i][k] * right[k][j];
					}
    			}
		}
		
		return result;
	}

	SquareMatrix<ClosestTwoDegree<Max<M, N, L>::value>::value, Field> A, B;

	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			A[i][j] = left[i][j];
		}
	}
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < L; ++j) {
			B[i][j] = right[i][j];
		}
	}

	A *= B;
	Matrix<M, L, Field> result;

	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < L; ++j) {
			result[i][j] = A[i][j];
		}
	}

	return result;
}

template <size_t M, size_t N, typename Field, size_t K, size_t L>
bool operator==(const Matrix<M, N, Field>& left, const Matrix<K, L, Field>& right) {
	if (M != K || N != L) {
		return false;
	}
	
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < M; ++j) {
			if (left[i][j] != right[i][j]) {
				return false;
			}
		}
	}
	
	return true;	
}
template <size_t M, size_t N, typename Field, size_t K, size_t L>
bool operator!=(const Matrix<M, N, Field>& left, const Matrix<K, L, Field>& right) {
	if (M != K || N != L) {
		return true;
	}
	
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < M; ++j) {
			if (left[i][j] != right[i][j]) {
				return true;
			}
		}
	}
	
	return false;
}
