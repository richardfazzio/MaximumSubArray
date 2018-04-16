# MaximumSubArray
This Code contains the solution for the Maximum Subarray problem in C++ both using BRUTE-FORCE and DIVIDE-AND-CONQUER.
The Maximum Subarray problem is finding the maximum value sum of a subarray of a given array as well as the indices of the subarray.
The practical applications of this is if you were analyzing changes over time, in say the stock market, and need to know when would 
be the optimal time to buy and sell. This algorithm will determine those points as well as the value of the given interval. BRUTE-FORCE
is a naive approach that simply goes through all possibilities and finds the max value RT=O(n^2), where n is the number of elements
of the original array. DIVIDE-AND-CONQUER breaks the problem into sub problems ( splitting each subarray into two, so on and so forth...)
then deciding which subarray contains a larger value, as well as the cross section of the subarray must be taken into account RT=O(nlgn).
It is also set up to do this experiment between the values of '-bounds-1' to '+bounds'.
It will perform each run 'm' times then take the average runtime of each 'run'.
The input sizes range from 5000,10000,...,100000 element arrays.
