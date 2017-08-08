all: erase compile

erase:
	rm -f a-helix c-helix d-helix
compile:
	cc -o a-helix all_helix.c -lm

	cc -o c-helix chain_helix.c -lm

	cc -o d-helix dom_helix.c -lm
