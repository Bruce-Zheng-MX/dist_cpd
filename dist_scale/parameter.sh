#!/bin/bash
for a in `seq 0 0.04 0.4`; 
  do 
    for b in {1..500};
      do
        echo "$a $b"
      done;
  done;
