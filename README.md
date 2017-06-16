### Dayhoff and the origin of PAM matrices

I was notified by email about a comment on an old [post](http://telliott99.blogspot.com/2008/08/pam-point-accepted-mutation.html)  (from 2008).  In the post, I discussed PAM matrices and Margaret Dayhoff.

As with all the old posts, the original linked files disappeared from Apple's servers long ago, but I was able to dig out archived material from an old disk drive to look for files related to Dayhoff.

I will try to reconstruct how it all works here.  Since the discussion is nearly ten years old, I decided not to try to fix a few errors that I came across in reconstructing the paper's tables.

#### Dayhoff et al. 1978

A big issue has always been poor quality for the scan of the paper that I found on the web.  The one I got through inter-library loan was no better (actually, it seems identical).  

There is a better scan available now:

http://compbio.berkeley.edu/class/c246/Reading/dayhoff-1978-apss.pdf

In particular, Fig 80 is very clear in the new version.  I put copies of these papers in the repo (``original papers`` directory), I  trust I won't get a takedown notice out of it.

#### original data

I found (many copies) of an old data file in my search of an archival hard drive.

* [dayhoff.corrected.changes.txt](data/dayhoff.corrected.changes.txt)

I don't recall the source of the corrections, but the file consists of a lower triangular matrix of observed changes between amino acids in the proteins catalogued by Dayhoff and colleagues.

There is also appended what appears to be the outcome of a ``diff`` comparison of the corrected and original versions.
 
I decided to reformat this file with right-justified values, and spaces rather than tabs.  I wrote a short Python [script](script.py) for this.  The output was written to [dayhoff.fig80.txt](data/dayhoff.fig80.txt).

I added a last row with the amino acids in order.

I noticed 2 big errors in comparison with the cleaned up pdf of the paper that I found on the web.  I read the original as KR=417 but the clear copy is 477.  I also read the original as VP=150 but the clear copy shows 50.

Here is a [screenshot](screenshots/Dayhoff Fig 80.png) of Fig 80 from the original paper.

### Redoing the analysis

#### Loading the data

The main Python script for working with the data is [dayhoff.py](dayhoff.py).  I'm going to work through the steps in the code and document them here.

The order of the amino acids is a little strange.  They are in lexicographical order for the 3-letter abbreviations, but the symbols I'm using are 1-letter.  Here is the order:

```
A R N D C Q E G H I L K M F P S T W Y V
```

The script saves that last row to [dayhoff.aaorder.txt](data/dayhoff.aaorder.txt), if the function ``save_order()`` is uncommented.

The first part of the code constructs a dictionary to hold the changes between amino acids in ``changeD``.  Changes are entered symmetrically, the value for **WD** is equal to the value for **DW**.

Next, the table is printed for proofreading.  Here is a [screenshot](screenshots/Screen Shot 2017-06-16 at 12.05.59 PM.png) of the output, which could be proofed against the input.

At this point we also tabulate the total number of changes observed involving the most observed changes (alanine, 3644 changes) compared to the least (tryptophan, 79 changes).  There is an important normalization that is needed before we can infer anything about mutability.

#### Origin of Fig 80

What Dayhoff et al did was to take 71 groups of closely related proteins known in 1978.  They observed a total of 1572 amino acid changes in these groups.  

They also inferred ancestral sequences:  i.e., given a tree (how made?) (ACGH, DBGH), (ADIJ, CBIJ)
infer ancestral sequences (ABGH, ABIJ).  The method for constructing the tree is not specified, but I believe they probably used parsimony.

They show this generic example:

```
ABGH -> ACGH and DBGH, infer A -> D and B -> C
ABIJ -> ADIJ and CBIJ, infer A -> C and B -> D
```

The first ancestor not given but changes are 

```
G <-> I and J <-> H.
```

From this analysis, they make a table where AC and CA both get 1, etc.  Note that CD is 0, even though same position (2) has both C and D, because it is inferred that this change was not direct, based on the tree.

The total of all the values in the dictionary is 31430 (next section).  Each change is counted twice, dividing by two we obtain 15715.  As described in the paper the total of all changes was 1572, the legend to Fig 80 explains that values were multiplied by 10.  The occurence of values in the table that are not multiples of 10 is explained by the comment

> (In practice, some of thepositions in the nodal sequencesare blank [ambiguous]. For thesew, e have treatedthe changes statistically,distributing them among all observed alternatives.)

#### Relative mutability

In order to say anything about the propensity of one amino acid to change into another (or simply to change at all), we need to know the frequencies of appearance of amino acids in the proteins Dayhoff et al studied.  This is provided in [Table 22](screenshots/Screen Shot 2017-06-16 at 12.26.44 PM.png).  I transcribed these as [dayhoff.frequencies.txt](data/dayhoff.frequencies.txt).  These were rewritten in row format for use with R ([dayhoff.frequencies.v2.txt](data/dayhoff.frequencies.v2.txt)).

Mutability is simply the ratio of observed changes to amino acid frequency, normalized to give ``Alanine = 100``.  Here is a [screenshot](screenshots/Screen Shot 2017-06-16 at 12.33.10 PM.png) of the output, which should be compared with [Table 21](screenshots/Dayhoff Table 21.png).

There are small discrepancies (compare 

* N at 134 v. 135
* S at 120 v. 119
* E at 102 v. 101
* T at 98 v. 97
* H at 66 v. 65
* W at 18 v. 19

I'm not sure of the reason for this, but we've shown we can get very close to their values.  From this point on, we use their values in [dayhoff.mutabilities.txt](data/dayhoff.mutabilities.txt).

#### PAM1 matrix

The last step of this first part of the analysis is to generate the file [PAM1.txt](data/PAM1.txt), which should be compared with [Fig 82](screenshots/Dayhoff Fig 82.png) of the text.

> An element of this matrix, Mij, gives the probability that the amino acid in column j will be replaced by the amino acide in row i after a given evolutionary interval, in this case 1 PAM
> 

As [described](screenshots/computation.png) in the paper, the off-diagonal entries in the PAM1 matrix are computed from the data in Fig 80.  

The original amino acids are specified in the columns (amino acid j), and the replacement amino acids are rows (amino acid i).

The first step in computation of the ij-th entry in PAM1 from the ij-th entry in Fig 80 is to form the ratio Aij/sum over i(Aij).  This is the ratio of observed changes involving amino acid i and j, divided by all changes involving amino acid j.

An additional factor is the mutability of amino acid j.  This seems kind of strange, since the mutability is sum over i(Aij) divided by the frequency of j, fj.  In other words, sum over i(Aij) cancels, leaving Aij/fj

The last factor is a proportionality constant, &lambda;, that is the same for all columns.

> The quantity 100 X ZfiMii gives the number of amino acids that will remain unchanged when a protein 100 links long, of average composition, is exposed to theevolu- tionary change represented by this matrix. This apparent evolutionary change depends upon thechoice of X, in this case ,chosen so that this change is 1 mutation. Since there are almost no superimposedchanges, this alsorepresents 1 PAM of change. If h hadbeen four times as large, the initial matrix would have represented 4 PAMs;  the discussion which follows would not be changed appreciably.

I didn't proof this carefully, but there are some differences, which are small, and I believe are contained to the diagonal.  The on-diagonal values should be 10000 minus the sum of the rest of that column.  The error is in my code.  I believe the problem is that I computed the total of the other entries before rounding the result.  They rounded the entries, then computed the total to subtract from 10000.

#### Continuing

You should now be able to follow the next [post](http://telliott99.blogspot.com/2008/08/pam-projecting-in-time.html) to generate the PAM250 matrix and so on.

If I get to it I will try to show how that part is done, but it is basically a matter of following the code in [Dayhoff.Rcode.txt](Dayhoff.Rcode.txt).