#------------------------------------------------------------------------------#
# Unit test configuration
#------------------------------------------------------------------------------#

bld.test('harris', 'harris.cxx ${top_srcdir}/src/main.cxx ${top_srcdir}/src/deck_wrapper.cxx')
