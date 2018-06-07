#!/bin/bash

# code formatter
for f in ./*.cpp ./*.h
do
	mv "$f" "$f"_
	/usr/bin/clang-format "$f"_ > "$f"
	rm -f "$f"_
done
