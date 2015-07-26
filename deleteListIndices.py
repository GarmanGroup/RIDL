# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 13:52:51 2015

@author: charlie
"""
from sys import stdout

# simple function to remove specific indices of list
def multi_delete(list_, ind2remove):
    ind2remove.sort(reverse=True)
      
    # unessential loading bar add-in (part1)
    listlen = len(ind2remove)
    point = listlen / 100
    increment = listlen / 100    
    
    counter = 0
    for index in ind2remove:
        del list_[index]
        
        # unessential loading bar add-in (part2)
        if(counter % (1 * point) == 0):
            stdout.write(
                    "\r[" + "=" * (counter / increment) +
                    " " * ((listlen - counter)/ increment) + 
                    "]" +  str(counter / point) + "%")
            stdout.flush()
        counter += 1
        
    return list_  