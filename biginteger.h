#include <iostream>
#include <vector>
#include <string>

class BigInteger {
  static const int64_t kBase = 1000000000;
  static void appendWithNulls(std::string &s);  // for output
  static int64_t makeBigDigit(std::string s);  // for input

 private:
  std::vector<int64_t> body_;
  int8_t sign_;  //  0 is null, 1 is positive, -1 if negative
  
  void killLeadingNulls();
  void shiftRight();  // for division

 public:
  BigInteger();
  BigInteger(int num);
  explicit BigInteger(std::string s);
  explicit operator bool() const;

  BigInteger operator+() const;
  BigInteger operator-() const;
  BigInteger &operator+=(const BigInteger &right);
  BigInteger &operator-=(const BigInteger &right);
  BigInteger &operator*=(const BigInteger &right);
  BigInteger &operator/=(const BigInteger &right);
  BigInteger &operator%=(const BigInteger &right);
  BigInteger &operator++();
  BigInteger &operator--();
  BigInteger operator++(int);
  BigInteger operator--(int);
  BigInteger gcd(BigInteger left, BigInteger right);
  BigInteger squareRoot(const BigInteger& num);

  std::string toString() const;
  bool isEven() const;
  void halve();
  bool isNegative() const;
  const std::vector<int64_t>& getBody() const;  // safe access to 'body vector'
  void changeSign();
};

bool operator==(const BigInteger &left, const BigInteger &right);
bool operator<(const BigInteger &left, const BigInteger &right);
bool operator!=(const BigInteger &left, const BigInteger &right);
bool operator<=(const BigInteger &left, const BigInteger &right);
bool operator>(const BigInteger &left, const BigInteger &right);
bool operator>=(const BigInteger &left, const BigInteger &right);

BigInteger operator+(const BigInteger &left, const BigInteger &right);
BigInteger operator-(const BigInteger &left, const BigInteger &right);
BigInteger operator*(const BigInteger &left, const BigInteger &right);
BigInteger operator/(const BigInteger &left, const BigInteger &right);
BigInteger operator%(const BigInteger &left, const BigInteger &right);

std::ostream& operator<<(std::ostream &out, const BigInteger &num);
std::istream& operator>>(std::istream &in, BigInteger &num);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void BigInteger::appendWithNulls(std::string& s) {
  size_t size = s.size();
  if (size < 9) {
    std::string temp = "0";
    for (unsigned int i = 1; i < 9 - size; ++i) {
      temp += '0';
    }
    s = temp + s;
  }
}
int64_t BigInteger::makeBigDigit(std::string s) {
  int64_t answer = 0;
  int64_t ten = 1;
  for (int i = s.size() - 1; i >= 0; --i) {
    answer += ten * (s[i] - '0');
    ten *= 10;
  }

  return answer;
}

void BigInteger::killLeadingNulls() {
  while (!body_.empty()) {
    if (body_.back() == 0) {
      body_.pop_back();
    } else break;
  }
  if (body_.empty()) {
    body_.push_back(0);
    sign_ = 0;
  }
}
void BigInteger::shiftRight() {
  if (sign_ == 0) {
    return;
  }
  
  body_.push_back(body_[body_.size() - 1]);
  for (int i = body_.size() - 2; i > 0; --i) {
    body_[i] = body_[i - 1];
  }
  body_[0] = 0;
}

BigInteger::BigInteger() : sign_(0) {
  body_.push_back(0);
}
BigInteger::BigInteger(int num) : sign_(1) {
  if (num == 0) {
    sign_ = 0;
    body_.push_back(0);
    return;
  } 
  if (num < 0) {
    num *= -1;
    sign_ = -1;
  }

  do {
    body_.push_back((num % kBase));
    num /= kBase;
  } while (num >= kBase);
  if (num != 0) {
    body_.push_back(num);
  }
}
BigInteger::BigInteger(std::string s) : sign_(1) {
  if (s[0] == '0') {
    sign_ = 0;
    body_.push_back(0);
    return;
  } 
  if (s[0] == '-') {
    sign_ = -1;
    s = s.substr(1);
  } else if (s[0] == '+') {
    s = s.substr(1);
  }

  for (int i = s.length(); i > 0; i -= 9) {
    if (i < 9) {
      body_.push_back(makeBigDigit(s.substr(0, i)));
    } else {
      body_.push_back(makeBigDigit(s.substr(i - 9, 9)));
    }
  }
  killLeadingNulls();
}
BigInteger::operator bool() const {
  return sign_ != 0;
}

BigInteger BigInteger::operator+() const {
  return {*this};
}
BigInteger BigInteger::operator-() const {
  BigInteger temp(*this);
  temp.changeSign();
  return temp;
}
BigInteger &BigInteger::operator+=(const BigInteger &right) {
  if (sign_ == 0) {
    *this = right;
    return *this;
  }
  if (right.sign_ == 0) {
    return *this;
  }
  if (sign_ < right.sign_) {
    BigInteger temp = right;
    temp -= -(*this);
    *this = temp;
    return *this;
  }
  if (sign_ > right.sign_) {
    *this -= -right;
    return *this;
  }

  int64_t carry = 0;
  for (unsigned int i = 0; i < std::max(body_.size(), right.body_.size()) || carry; ++i) {
    if (i == body_.size()) {
      body_.push_back(0);
    }
    body_[i] += carry + (i < right.body_.size() ? right.body_[i] : 0);
    carry = body_[i] >= kBase;
    if (carry != 0) {
      body_[i] -= kBase;
    }
  }

  return *this;
}
BigInteger &BigInteger::operator-=(const BigInteger &right) {
  if (sign_ == 0) {
    (*this) = -right;
    return *this;
  }
  if (right.sign_ == 0) {
    return *this;
  }
  if (sign_ != right.sign_) {
    *this += -right;
    return *this;
  }
  if (sign_ == -1) {
    BigInteger pos_this = -(*this);
    if (-right >= pos_this) {
      BigInteger temp = -right;
      temp -= pos_this;
      *this = temp;
      return *this;
    }
    pos_this -= -right;
    *this = -pos_this;
    return *this;
  }
  if (sign_ == 1 && *this < right) {
    BigInteger temp = right;
    temp -= *this;
    *this = -temp;
    return *this;
  }

  int64_t carry = 0;
  for (unsigned int i = 0; i < body_.size() || carry; ++i) {
    body_[i] -= carry + (i < right.body_.size() ? right.body_[i] : 0);
    carry = body_[i] < 0;
    if (carry != 0) {
      body_[i] += kBase;
    }
  }
  killLeadingNulls();

  return *this;
}
BigInteger &BigInteger::operator*=(const BigInteger &right) {
  if (sign_ == 0) {
    return *this;
  }
  if (right.sign_ == 0) {
    sign_ = 0;
    body_ = right.body_;
    return *this;
  }

  std::vector<int64_t> answer(body_.size() + right.body_.size());
  int64_t carry;
  for (unsigned int i = 0; i < body_.size(); ++i) {
    carry = 0;
    for (unsigned int j = 0; j < right.body_.size() || carry != 0; ++j) {
      int64_t cur = answer[i + j] + body_[i] * (j < right.body_.size() ? right.body_[j] : 0) + carry;
      answer[i + j] = cur % kBase;
      carry = cur / kBase;
    }
  }
  body_ = answer;
  sign_ *= right.sign_;
  killLeadingNulls();

  return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &right) {
  if (sign_ == 0) {
    return *this;
  }
  
  if (right.isNegative()) {
    changeSign();
    *this /= -right;
    return *this;
  }

  if (right == 2) {
    halve();
    return *this;
  }

  BigInteger res;
  res.body_.resize(body_.size());
  BigInteger curValue;

  for (int i = body_.size() - 1; i >= 0; --i) {
    curValue.shiftRight();
    curValue += body_[i];
    int64_t x = 0;
    int64_t l = 0, r = kBase;
    while (l <= r) {
      int64_t m = (l + r) >> 1;
      BigInteger cur = right * m;
      if (cur <= curValue) {
        x = m;
        l = m + 1;
      } else {
        r = m - 1;
      }
    }

    res.body_[i] = x;
    curValue = curValue - right * x;
  }

  sign_ *= right.sign_;
  body_ = res.body_;
  killLeadingNulls();

  return *this;
}
BigInteger &BigInteger::operator%=(const BigInteger &right) {
  if (sign_ == 0) {
    return *this;
  }

  if (right.isNegative()) {
    changeSign();
    *this %= -right;
    return *this;
  }

  BigInteger res;
  res.body_.resize(body_.size());
  BigInteger curValue;

  for (int i = body_.size() - 1; i >= 0; --i) {
    curValue.shiftRight();
    curValue += body_[i];
    int64_t x = 0;
    int64_t l = 0, r = kBase;
    while (l <= r) {
      int64_t m = (l + r) >> 1;
      BigInteger cur = right * m;
      if (cur <= curValue) {
        x = m;
        l = m + 1;
      } else {
        r = m - 1;
      }
    }

    res.body_[i] = x;
    curValue = curValue - right * x;
  }

  body_ = curValue.body_;
  killLeadingNulls();

  return *this;
}
BigInteger &BigInteger::operator++() {
  *this += 1;
  return *this;
}
BigInteger &BigInteger::operator--() {
  *this -= 1;
  return *this;
}
BigInteger BigInteger::operator++(int) {
  BigInteger temp = *this;
  ++(*this);
  return temp;
}
BigInteger BigInteger::operator--(int) {
  BigInteger temp = *this;
  --(*this);
  return temp;
}
BigInteger Gcd(BigInteger left, BigInteger right) {
  if (left.isNegative()) {
    left.changeSign();
  }
  if (right.isNegative()) {
    right.changeSign();
  }
  BigInteger gcd = 1;
  while (left != 0 && right != 0) {
    if (left.isEven() && right.isEven()) {
      left.halve();
      right.halve();
      gcd += gcd;
    } else if (left.isEven()) {
      left.halve();
    } else if (right.isEven()) {
      right.halve();
    } else {
      if (left < right) {
        BigInteger temp = right;
        right = left;
        left = temp;
      }
      left -= right;
    }
  }

  return (left == 0 ? right : left) * gcd;
}
BigInteger SquareRoot(const BigInteger &num) {	
  if (num == 0 || num == 1) {
    return num;
  }

  BigInteger r = num;
  BigInteger l = 1;
  BigInteger m, temp;

  while (l < r - 1) {
    m = l + r;
	m.halve();
    temp = m * m;

    if (num == temp) {
		return m;
	}
	if (num < temp) {
      r = m;
    } else {
      l = m;
	}
  }

  return l;
}

std::string BigInteger::toString() const {
  std::string s = std::to_string(body_[body_.size() - 1]);
  if (sign_ == -1) {
    s = '-' + s;
  }
  if (body_.size() == 1) {
    return s;
  }

  std::string temp;
  for (int i = body_.size() - 2; i >= 0; --i) {
    temp = std::to_string(body_[i]);
    appendWithNulls(temp);
    s += temp;
  }

  return s;
}
bool BigInteger::isEven() const {
  return (body_[0] % 2 == 0);
}
void BigInteger::halve() {
  int64_t carry = 0;
  int64_t current;
  for (int i = body_.size() - 1; i >= 0; --i) {
    current = body_[i] + carry * kBase;
    body_[i] = current / 2;
    carry = current % 2;
  }

  killLeadingNulls();
}
bool BigInteger::isNegative() const {
  return sign_ == -1;
}
const std::vector<int64_t>& BigInteger::getBody() const {
  return body_;
}
void BigInteger::changeSign() {
    sign_ *= -1;
}

BigInteger operator+(const BigInteger &left, const BigInteger &right) {
  BigInteger tmp = left;
  tmp += right;
  return tmp;
}
BigInteger operator-(const BigInteger &left, const BigInteger &right) {
  BigInteger tmp = left;
  tmp -= right;
  return tmp;
}
BigInteger operator*(const BigInteger &left, const BigInteger &right) {
  BigInteger tmp = left;
  tmp *= right;
  return tmp;
}
BigInteger operator/(const BigInteger &left, const BigInteger &right) {
  BigInteger tmp = left;
  tmp /= right;
  return tmp;
}
BigInteger operator%(const BigInteger &left, const BigInteger &right) {
  BigInteger tmp = left;
  tmp %= right;
  return tmp;
}

bool operator==(const BigInteger &left, const BigInteger &right) {
  return left.getBody() == right.getBody() && left.isNegative() == right.isNegative();
}
bool operator<(const BigInteger &left, const BigInteger &right) {
  if (left == right) {
    return false;
  }
  if (left.isNegative()) {
    if (right.isNegative()) {
      return ((-right) < (-left));
    }
    return true;
  }
  if (right.isNegative()) {
    return false;
  }
  if (left.getBody().size() != right.getBody().size()) {
    return left.getBody().size() < right.getBody().size();
  }

  for (int i = left.getBody().size() - 1; i >= 0; --i) {
    if (left.getBody()[i] != right.getBody()[i]) {
      return left.getBody()[i] < right.getBody()[i];
    }
  }

  return false;
}
bool operator!=(const BigInteger &left, const BigInteger &right) {
  return !(left == right);
}
bool operator<=(const BigInteger &left, const BigInteger &right) {
  return (left < right || left == right);
}
bool operator>(const BigInteger &left, const BigInteger &right) {
  return !(left <= right);
}
bool operator>=(const BigInteger &left, const BigInteger &right) {
  return !(left < right);
}

std::ostream &operator<<(std::ostream &out, const BigInteger &num) {
  out << num.toString();
  return out;
}
std::istream &operator>>(std::istream &in, BigInteger &num) {
  std::string s;
  in >> s;
  num = BigInteger{s};
  return in;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Rational {
 private:
  BigInteger num_;
  BigInteger den_;

  void toNormal();

 public:
  Rational();
  Rational(int n);
  Rational(BigInteger n);
  Rational(BigInteger n, BigInteger d);
  Rational(const Rational &num);
  Rational &operator=(const Rational &num);
  ~Rational() = default;
  explicit operator double() const;

  const BigInteger &getNumerator() const;
  const BigInteger &getDenominator() const;
  std::string toString() const;
  std::string asDecimal(size_t precision = 0) const;

  Rational &operator+=(const Rational &right);
  Rational &operator-=(const Rational &right);
  Rational &operator*=(const Rational &right);
  Rational &operator/=(const Rational &right);
  Rational &operator++();
  Rational &operator--();
  Rational operator++(int);
  Rational operator--(int);
  Rational operator-() const;
};

Rational operator+(const Rational &right, const Rational &left);
Rational operator-(const Rational &right, const Rational &left);
Rational operator*(const Rational &right, const Rational &left);
Rational operator/(const Rational &right, const Rational &left);
bool operator<(const Rational &right, const Rational &left);
bool operator>(const Rational &right, const Rational &left);
bool operator<=(const Rational &right, const Rational &left);
bool operator>=(const Rational &right, const Rational &left);
bool operator==(const Rational &right, const Rational &left);
bool operator!=(const Rational &right, const Rational &left);

std::ostream &operator<<(std::ostream &out, const Rational &r);
std::istream &operator>>(std::istream &in, Rational &r);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Rational::toNormal() {
  BigInteger gcd = Gcd(num_, den_);
  den_ /= gcd;
  num_ /= gcd;
  if (den_ < 0) {
    den_.changeSign();
    num_.changeSign();
  }
}

Rational::Rational() : num_(0), den_(1) {
}
Rational::Rational(int n) : num_(n), den_(1) {
}
Rational::Rational(BigInteger n) : num_(n), den_(1) {
}
Rational::Rational(BigInteger n, BigInteger d) : num_(n), den_(d) {
  toNormal();
}
Rational::Rational(const Rational &num) {
  num_ = num.num_;
  den_ = num.den_;
  toNormal();
}
Rational &Rational::operator=(const Rational &num) {
  if (&num == this) {
    return *this;
  }

  num_ = num.num_;
  den_ = num.den_;
  toNormal();

  return *this;
}
Rational::operator double() const {
  return std::stod(asDecimal(16));
}

const BigInteger &Rational::getNumerator() const {
  return num_;
}
const BigInteger &Rational::getDenominator() const {
  return den_;
}
std::string Rational::toString() const {
  std::string s = num_.toString();
  if (den_ != 1) {
    s += '/' + den_.toString();
  }

  return s;
}
std::string Rational::asDecimal(size_t precision) const {
  if (precision == 0) {
    return (num_ / den_).toString();
  }

  if (den_ == 1) {
    return (num_ / den_).toString() + '.' + std::string(precision, '0');
  }

  std::string s = '1' + std::string(precision + 1, '0');
  BigInteger temp = num_;
  temp *= BigInteger(s);
  temp /= den_;

  if (temp.getBody()[0] % 10 >= 5) {
    temp += (temp.isNegative() ? -10 : 10);
  }

  if (temp.isNegative()) {
    s = (-temp).toString();
  } else {
    s = temp.toString();
  }
  s.pop_back();

  if (s.length() <= precision) {
    return (temp.isNegative() ? "-0." : "0.") + std::string(precision - s.length(), '0') + s;
  }
  if (temp.isNegative()) {
    return '-' + s.substr(0, s.length() - precision) + '.' + s.substr(s.length() - precision);
  }
  return s.substr(0, s.length() - precision) + '.' + s.substr(s.length() - precision);
}

Rational &Rational::operator+=(const Rational &right) {
  num_ *= right.den_;
  num_ += right.num_ * den_;
  den_ *= right.getDenominator();
  toNormal();
  return *this;
}
Rational &Rational::operator-=(const Rational &right) {
  num_ *= right.getDenominator();
  num_ -= right.getNumerator() * den_;
  den_ *= right.getDenominator();
  toNormal();
  return *this;
}
Rational &Rational::operator*=(const Rational &right) {
  num_ *= right.getNumerator();
  den_ *= right.getDenominator();
  toNormal();
  return *this;
}
Rational &Rational::operator/=(const Rational &right) {
  num_ *= right.getDenominator();
  den_ *= right.getNumerator();
  toNormal();
  return *this;
}
Rational &Rational::operator++() {
  num_ += den_;
  return *this;
}
Rational &Rational::operator--() {
  num_ -= den_;
  return *this;
}
Rational Rational::operator++(int) {
  Rational temp(num_, den_);
  ++(*this);
  return temp;
}
Rational Rational::operator--(int) {
  Rational temp(num_, den_);
  --(*this);
  return temp;
}
Rational Rational::operator-() const {
  return {-num_, den_};
}
Rational operator+(const Rational &right, const Rational &left) {
  BigInteger n = right.getNumerator() * left.getDenominator() + left.getNumerator() * right.getDenominator();
  BigInteger d = right.getDenominator() * left.getDenominator();
  return {n, d};
}
Rational operator-(const Rational &right, const Rational &left) {
  BigInteger n = right.getNumerator() * left.getDenominator() - left.getNumerator() * right.getDenominator();
  BigInteger d = right.getDenominator() * left.getDenominator();
  return {n, d};
}
Rational operator*(const Rational &right, const Rational &left) {
  BigInteger n = right.getNumerator() * left.getNumerator();
  BigInteger d = right.getDenominator() * left.getDenominator();
  return {n, d};
}
Rational operator/(const Rational &right, const Rational &left) {
  BigInteger n = right.getNumerator() * left.getDenominator();
  BigInteger d = right.getDenominator() * left.getNumerator();
  return {n, d};
}
bool operator<(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r < l);
}
bool operator>(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r > l);
}
bool operator<=(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r <= l);
}
bool operator>=(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r >= l);
}
bool operator==(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r == l);
}
bool operator!=(const Rational &right, const Rational &left) {
  BigInteger r = right.getNumerator() * left.getDenominator();
  BigInteger l = right.getDenominator() * left.getNumerator();
  return (r != l);
}

std::ostream &operator<<(std::ostream &out, const Rational &r) {
  out << r.toString();
  return out;
}

std::istream &operator>>(std::istream &in, Rational &r) {
  std::string s;
  bool flag = true;
  in >> s;

  for (unsigned int i = 0; i < s.length(); ++i) {
    if (s[i] == '/') {
      std::string s_num, s_den;
      s_num = s.substr(0, i);
      s_den = s.substr(i + 1, s.length() - i - 1);
      r = Rational(BigInteger(s_num), BigInteger(s_den));
      flag = false;
      break;
    }
  }

  if (flag) {
    r = Rational(BigInteger(s), 1);
  }

  return in;
}
