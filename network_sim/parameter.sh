#!/bin/bash
for a in `seq 0.0000000001 0.05 0.6000000001`; 
  do 
    for b in {1..500};
      do
        echo "$a $b"
      done;
  done;
