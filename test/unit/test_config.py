#------------------------------------------------------------------------------#
# Unit test configuration
#------------------------------------------------------------------------------#

bld.test('accel',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=accel.deck', [1])
bld.test('cyclo',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=cyclo.deck', [1])
bld.test('inbndj',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=inbndj.deck', [1])
bld.test('interpe',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=interpe.deck', [1])
bld.test('outbndj',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=outbndj.deck', [1])
bld.test('pcomm',
	'${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx',
	'-DINPUT_DECK=pcomm.deck', [8])
