// Macro that is to be included inside advance_p pipelined variants.
// In an ideal world this would not exist, and the pipeline code would not be
// repeated in multiple files
// This is done so this common element can appear only once, instead of once
// per file
#define MOVE_OUTBND(N)                                              \
if ( outbnd(N) )                                /* Unlikely */      \
{                                                                   \
    local_pm->dispx = ux(N);                                          \
    local_pm->dispy = uy(N);                                          \
    local_pm->dispz = uz(N);                                          \
    local_pm->i     = ( p - p0 ) + N;                                 \
    if ( move_p( p0, local_pm, a0, g, _qsp ) )    /* Unlikely */      \
    {                                                                 \
        if ( nm < max_nm )                                              \
        {                                                               \
            copy_4x1( &pm[nm++], local_pm );                              \
        }                                                               \
        else                                        /* Unlikely */      \
        {                                                               \
            itmp++;                                                       \
          /* Also undo the shift that move_p did, to keep p->i in a valid range */ \
          /* If we got here, we're running the risk of ruining the physics of  */\
          /* the simulation. Take the mover warning very seriously. */ \
          p0[ local_pm->i ].i = p0[ local_pm->i].i >> 3; \
        }                                                               \
    }                                                                 \
}
