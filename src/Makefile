ifndef ACGTUN_UTIL
ACGTUN_UTIL=$(AFFINE-GAP-PENALTY-ALIGNMENT)/src/acgtun_util/
endif

app_subdirs=gap_alignment

all:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(AFFINE-GAP-PENALTY-ALIGNMENT) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(AFFINE-GAP-PENALTY-ALIGNMENT) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(AFFINE-GAP-PENALTY-ALIGNMENT) test; \
	done;

clean:
	@make -C acgtun_util OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} ACGTUN_UTIL=$(ACGTUN_UTIL) SRC_ROOT=$(AFFINE-GAP-PENALTY-ALIGNMENT) clean; \
	done;
