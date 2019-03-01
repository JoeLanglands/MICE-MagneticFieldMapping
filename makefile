SRCDIR = $(addprefix $(PWD), /micemag/makefields/)

default:
	$(MAKE) -C $(SRCDIR)


clean:
	$(MAKE) -C $(SRCDIR) clean

.PHONY: default clean
