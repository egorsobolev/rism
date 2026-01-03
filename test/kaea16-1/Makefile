test:	rismaw

mdp_10_amb.mdcrd:
	gzip -d -k mdp_10_amb.mdcrd.gz

kaea1_straj.bin: mdp_10_amb.mdcrd
	../../bin/straj -a 231 mdp_10_amb.mdcrd kaea1_straj.bin

kaea1_hist.bin:	kaea1_straj.bin
	../../bin/mkhist kaea1_straj.bin kaea1_hist.bin

aw-05-2048-056.bin:	kaea1_hist.bin
	../../bin/mkavgw kaea1_hist.bin -s 0.05:2048:0.056 kaea1_aw

rismaw:	aw-05-2048-056.bin
	../../bin/awrism -p kaea1.rtxt \
	                 -g kaea1_aw-05-2048-056.bin:../tip3p-dr056-np2048.rbin \
	                 kaea1

.PHONY: clean
clean:
	rm -f kaea1_aw-05-2048-056.bin
	rm -f kaea1_hist.bin
	rm -f kaea1_straj.bin
	rm -f mdp_10_amb.mdcrd
	rm -f tuv.bin
