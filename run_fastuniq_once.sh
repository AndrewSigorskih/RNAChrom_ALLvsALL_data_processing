#!/bin/bash
{ read LEFT; read RIGHT; } < list.tmp
fastuniq -i list2 -o fastuniq/${LEFT} -p fastuniq/${RIGHT}
echo "finished"
