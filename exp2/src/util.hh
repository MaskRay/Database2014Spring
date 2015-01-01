class CRC32Hash {
  size_t operator()(uint64_t x) const {
    return _mm_crc32_u64(0,x);
  }
  size_t operator()(int64_t x) const {
    return _mm_crc32_u64(0,x);
  }
  size_t operator()(uint32_t x) const {
    return _mm_crc32_u32(0,x);
  }
  size_t operator()(int32_t x) const {
    return _mm_crc32_u32(0,x);
  }

  inline unsigned operator()(const std::string& s) const {
    unsigned l=s.length();
    const char* from = s.data();
    if (l<=8) {
      return _mm_crc32_u64(0,(reinterpret_cast<const uint64_t*>(from)[0])<<(64-8*l));
    } else if (l<=16) {
      return _mm_crc32_u64(_mm_crc32_u64(0,reinterpret_cast<const uint64_t*>(from)[0]),(reinterpret_cast<const  uint64_t*>(from)[1])<<(128-8*l));
    } else if (l<=24) {
      return _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(0,reinterpret_cast<const uint64_t*>(from)[0]),reinterpret_cast<const  uint64_t*>(from)[1]),(reinterpret_cast<const uint64_t*>(from)[2])<<(128+64-8*l));
    } else {
      return _mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(_mm_crc32_u64(0,reinterpret_cast<const uint64_t*>(from)[0]),reinterpret_cast<const uint64_t*>(from)[1]),reinterpret_cast<const uint64_t*>(from)[2]),(reinterpret_cast<const uint64_t*>(from)[3])<<(256-8*l));
    }
  }
}
