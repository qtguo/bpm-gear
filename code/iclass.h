#ifndef __IHEAP_H__
#define __IHEAP_H__

const int VectorDefaultSize = 8;
template <typename _Key, typename _Value>
struct Key_Value {
  _Key key;
  _Value value;
  Key_Value(const _Key& k, const _Value& v) {
    key = k;
    value = v;
  }
  bool operator==(const Key_Value<_Key, _Value>& p) const { return (key == p.key); }
  bool operator<(const Key_Value<_Key, _Value>& p) const {
    if (key < p.key)
      return true;
    else if (key > p.key)
      return false;
    return (value < p.value);
  }

  Key_Value(int tmp) {
    key = tmp;
    value = -1;
  }

  Key_Value() {}
};

template <typename _T>
class iVector {
 public:
  int64_t m_size;
  int64_t m_num;
  int64_t m_valid_num;
  _T* m_data;

  void free_mem() { delete[] m_data; }

  iVector() {
    m_size = VectorDefaultSize;
    m_data = new _T[VectorDefaultSize];
    m_num = 0;
    m_valid_num = 0;
  }

  ~iVector() { free_mem(); }

  iVector(size_t n) {
    if (n == 0) {
      n = VectorDefaultSize;
    }
    //		printf("iVector allocate: %d\n",n);
    m_size = n;
    m_data = new _T[m_size];
    m_num = 0;
  }

  void push_back(const _T& d) {
    m_data[m_num] = d;
    m_num++;
    if (m_num == m_size) {
      re_allocate(m_size * 2);
    }
  }

  void* get_next_position() { return (void*)(m_data + m_num); }

  void increase_num() {
    m_num++;
    if (m_num == m_size) {
      re_allocate(m_size * 2);
    }
  }

  void push_back(const _T* p, int64_t len) {
    int64_t required_size = m_num + len + 1;
    while (required_size > m_size) {
      re_allocate(m_size * 2);
    }

    // memcpy(m_data + m_num, p, sizeof(_T) * len);
    std::copy(p, p + len, m_data + m_num);
    m_num += len;
  }

  void re_allocate(int64_t size) {
    if (size < m_num) {
      return;
    }
    _T* tmp = new _T[size];
    // memcpy(tmp, m_data, sizeof(_T) * m_num);
    std::copy(m_data, m_data + m_num, tmp);
    m_size = size;
    delete[] m_data;
    m_data = tmp;
  }

  void clean() {
    m_num = 0;
    m_valid_num = 0;
  }
  void assign(iVector& t) {
    m_num = t.m_num;
    m_size = t.m_size;
    delete[] m_data;
    m_data = t.m_data;
  }

  bool remove(_T& x) {
    for (int64_t l = 0, r = m_num; l < r;) {
      int64_t m = (l + r) / 2;

      if (m_data[m] == x) {
        m_num--;
        if (m_num > m) memmove(m_data + m, m_data + m + 1, sizeof(_T) * (m_num - m));
        return true;
      } else if (m_data[m] < x)
        l = m + 1;
      else
        r = m;
    }
    return false;
  }

  _T& operator[](int64_t i) { return m_data[i]; }
};

typedef iVector<rrid_t> InvList_t;

struct bArray {
  nid_t size;
  nid_t* arr;

  bArray() {
    size = 0;
    arr = 0;
  }

  void makeSpace(nid_t num) {
    size = (num >> 5) + 1;
    arr = new nid_t[size];
    for (nid_t i = 0; i < size; i++) {
      arr[i] = 0;
    }
  }

  inline nid_t* getData(nid_t index) { return arr + (index >> 5); }
  inline bool getBool(nid_t index) { return (arr[index >> 5]) & (1 << (index & 0x1F)); }
  inline void setTrue(nid_t index) { arr[index >> 5] |= (1 << (index & 0x1F)); }
  inline void setFalse(nid_t index) {
    nid_t flag = ~(1 << (index & 0x1F));
    arr[index >> 5] &= flag;
  }
};

#endif