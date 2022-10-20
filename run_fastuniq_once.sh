#!/bin/bash
{ read LEFT; read RIGHT; } < list.tmp
fastuniq -i list.tmp -o fastuniq/${LEFT} -p fastuniq/${RIGHT}
echo "finished"
