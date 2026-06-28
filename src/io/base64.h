#pragma once



#include <iosfwd>

template <typename T>
class Base64Encoder {
 private:
  std::ostream& ostr;
  size_t charFullLength;
  size_t numWritten;
  int numOverflow;
  unsigned char overflow[3];

  static constexpr char enc64[65] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz0123456789+/";

  void fillOverflow(const unsigned char* charData, size_t charLength, size_t& pos) {
    while (numOverflow < 3 && pos < charLength) {
      overflow[numOverflow] = charData[pos];
      ++numOverflow;
      ++pos;
    }
    if (numOverflow == 3) {
      encodeBlock(overflow);
      numOverflow = 0;
    }
  }
  void flushOverflow() {
    if (numOverflow > 0) {
      for (int iOverflow = numOverflow; iOverflow < 3; ++iOverflow) {
        overflow[iOverflow] = 0;
      }
      encodeUnfinishedBlock(overflow, numOverflow);
      numOverflow = 0;
    }
  }
  void encodeBlock(const unsigned char* data) {
    ostr << enc64[data[0] >> 2];
    ostr << enc64[((data[0] & 0x03) << 4) | ((data[1] & 0xf0) >> 4)];
    ostr << (unsigned char)(enc64[((data[1] & 0x0f) << 2) | ((data[2] & 0xc0) >> 6)]);
    ostr << (unsigned char)(enc64[data[2] & 0x3f]);
  }
  void encodeUnfinishedBlock(const unsigned char* data, int length) {
    ostr << enc64[data[0] >> 2];
    ostr << enc64[((data[0] & 0x03) << 4) | ((data[1] & 0xf0) >> 4)];
    ostr << (unsigned char)(length == 2 ? enc64[((data[1] & 0x0f) << 2) | ((data[2] & 0xc0) >> 6)]
                                        : '=');
    ostr << (unsigned char)('=');
  }

 public:
  Base64Encoder(std::ostream& ostr_, size_t fullLength_)
      : ostr(ostr_), charFullLength(fullLength_ * sizeof(T)), numWritten(0), numOverflow(0) {}

  void encode(const T* data, size_t length) {
    const unsigned char* charData = reinterpret_cast<const unsigned char*>(data);
    size_t charLength = length * sizeof(T);
    // numWritten+charLength <= charFullLength

    size_t pos = 0;
    fillOverflow(charData, charLength, pos);
    while (pos + 3 <= charLength) {
      encodeBlock(charData + pos);
      pos += 3;
    }
    fillOverflow(charData, charLength, pos);
    numWritten += charLength;
    if (numWritten == charFullLength) flushOverflow();
  }
};


template <typename T>
class Base64Decoder {
 private:
  void flushOverflow(unsigned char* charData, size_t charLength, size_t& pos) {
    while (posOverflow < 3 && pos < charLength) {
      charData[pos] = overflow[posOverflow];
      ++pos;
      ++posOverflow;
    }
  }
  unsigned char getNext() {
    unsigned char nextChar;
    istr >> nextChar;
    return (unsigned char)(dec64[nextChar - 43] - 62);
  }
  void decodeBlock(unsigned char* data) {
    unsigned char input[4];
    input[0] = getNext();
    input[1] = getNext();
    input[2] = getNext();
    input[3] = getNext();
    data[0] = (unsigned char)(input[0] << 2 | input[1] >> 4);
    data[1] = (unsigned char)(input[1] << 4 | input[2] >> 2);
    data[2] = (unsigned char)(((input[2] << 6) & 0xc0) | input[3]);
  }

  static constexpr char dec64[82] =
    "|###}rstuvwxyz{#######>?@"
    "ABCDEFGHIJKLMNOPQRSTUVW######XYZ"
    "[\\]^_`abcdefghijklmnopq";

  std::istream& istr;
  size_t charFullLength;
  size_t numRead;
  int posOverflow;
  unsigned char overflow[3];

 public:
  Base64Decoder(std::istream& istr_, size_t fullLength_)
      : istr(istr_), charFullLength(fullLength_ * sizeof(T)), numRead(0), posOverflow(3) {}
  void decode(T* data, size_t length) {
    unsigned char* charData = reinterpret_cast<unsigned char*>(data);
    size_t charLength = length * sizeof(T);
    // numRead+charLength <= charFullLength

    size_t pos = 0;
    flushOverflow(charData, charLength, pos);
    while (pos + 3 <= charLength) {
      decodeBlock(charData + pos);
      pos += 3;
    }
    if (pos < charLength) {
      decodeBlock(overflow);
      posOverflow = 0;
      flushOverflow(charData, charLength, pos);
    }
    numRead += charLength;
  }
};