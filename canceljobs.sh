#!/bin/bash

for id in `seq 38498454 1 38498483`;do
    scancel $id
done

