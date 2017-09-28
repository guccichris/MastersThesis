# This is the Lennard-Jones 12-6 potential equation for easy use in other code

# initialize the params w and sigma
w, sigma = 1, 1

# function to be imported like: from lennardJones12-6 import LJ
def LJ(r):
    return w*((sigma**12/r**12) - 2*(sigma**6/r**6))
