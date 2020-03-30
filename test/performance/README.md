# General Idea

The idea behind these tests is to give a user a quick way to benchmark/verify
their install of VPIC for a given platform, and to allow them to collect some
performance data in a reproducible and automated way.

# Tests

## Implemented

- Repeatedly center (`center_p`) and uncenter (`uncenter_p`) the particle
population to give an upper bound of "particle push" speed

## Future

- Add a test to do a full `advance_p` call and then undo it (`uncenter_p`?)

# TODO 

- Need to add print output to describe the platform (and possibly the install
        configuration of VPIC)

