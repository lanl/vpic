#ifndef BitField_h
#define BitField_h

class BitField
	{
	public:

		BitField(uint32_t setbits = 0xffffffff) : bits_(setbits) {}
		BitField(const BitField & bf) : bits_(bf.bits_) {}
		~BitField() {}

		/*!---------------------------------------------------------------------
		 * Set bits in mask.
		----------------------------------------------------------------------*/
		uint32_t set(uint32_t mask) {
			return bits_ |= mask;
		} // set

		/*!---------------------------------------------------------------------
		 * Set individual bit.
		----------------------------------------------------------------------*/
		uint32_t setbit(size_t bit) {
			uint32_t tmp = 1<<bit;
			return bits_ |= tmp;
		} // setbit

		/*!---------------------------------------------------------------------
		 * Clear bits in mask.
		----------------------------------------------------------------------*/
		uint32_t clear(uint32_t mask) {
			return bits_ &= ~mask;
		} // clear

		/*!---------------------------------------------------------------------
		 * Clear individual bit.
		----------------------------------------------------------------------*/
		uint32_t clearbit(size_t bit) {
			uint32_t tmp = 1<<bit;
			return bits_ &= ~tmp;
		} // setbit

		/*!---------------------------------------------------------------------
		 * Check whether bit at index is set.
		----------------------------------------------------------------------*/
		bool bitset(size_t bit) const {
			uint32_t tmp = 1<<bit;
			return tmp & bits_;
		} // setbit

		/*!---------------------------------------------------------------------
		 * Check whether bits in mask are set.
		----------------------------------------------------------------------*/
		bool bitsset(uint32_t mask) const {
			uint32_t tmp = 1<<mask;
			return tmp & bits_;
		} // setbit

		/*!---------------------------------------------------------------------
		 * Check whether bit at index is clear.
		----------------------------------------------------------------------*/
		bool bitclear(size_t bit) const {
			uint32_t tmp = 1<<bit;
			return !(tmp & bits_);
		} // bitclear

		/*!---------------------------------------------------------------------
		 * Return the sum of the set bits in the mask.
		----------------------------------------------------------------------*/
		size_t bitsum(const size_t * indeces, size_t size) const {
			size_t sum(0);
			for(size_t i(0); i<size; i++) {
				if(bitset(indeces[i])) {
					++sum;
				} // if
			} // for
			return sum;
		} // bitclear

		/*!---------------------------------------------------------------------
		 * Return the sum of the set bits.
		----------------------------------------------------------------------*/
		size_t bitsum() const {
			size_t sum(0);
			for(size_t i(0); i<32; i++) {
				if(bitset(i)) {
					++sum;
				} // if
			} // for
			return sum;
		} // bitclear

	private:

		uint32_t bits_;

	}; // class BitField

#endif // BitField_h
