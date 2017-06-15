### Dayhoff and the PAM matrices

I was notified by email about a comment on an old [post](http://telliott99.blogspot.com/2008/08/pam-point-accepted-mutation.html) on PAM matrices and Margaret Dayhoff (from 2008) along with a request for my code.

These files disappeared from Apple's servers long ago, but I was able to dig out archived material from an old disk drive to look for files related to Dayhoff.  Luckily the disk drive still works.

I have tried to reconstruct how it all worked.  It's challenging!

A big issue back then was the poor quality of the scan of the paper that I found on the web, and the one I got through inter-library loan was no better (actually, identical).  

There is a better scan available now:

http://compbio.berkeley.edu/class/c246/Reading/dayhoff-1978-apss.pdf

Fig 80 is very clear in the new version.  

So the first thing was to proof my file of these "accepted" changes against the new text.  Actually the very first thing is to reformat my file with right-justified values and spaces rather than tabs.  I notice 2 big errors in the original file.

I read the original as KR=417 but the clear copy is 477.  I also read the original as VP=150 but the clear copy shows 50.

This repository has the original data files as well as the reworked ones.  ``script.py`` is just to reformat the data in ``dayhoff.corrected.changes.txt``.  The new version is ``dayhoff.fig80.txt``.

The main Python script for working with the data is ``dayhoff.py``.  I hope it is self-explanatory.  There are intermediate steps (commented out) to print out versions of the data along the way.

The end result is ``PAM1.txt``, which is Fig 82 of the text.

You should be able to follow the next [post](http://telliott99.blogspot.com/2008/08/pam-projecting-in-time.html) to generate the PAM250 matrix and so on.

I hope this helps.