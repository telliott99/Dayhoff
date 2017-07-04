### Dayhoff and the origin of PAM matrices

I was notified by email about a comment on an old [post](http://telliott99.blogspot.com/2008/08/pam-point-accepted-mutation.html) from 2008.  There, I discussed PAM matrices and Margaret Dayhoff.

As with all the old posts, the original linked files disappeared from Apple's servers long ago, but I was able to dig out archived material from an old disk drive to look for files related to Dayhoff.

I will try to reconstruct how it all works here.  Since the discussion is nearly ten years old, I decided not to try to troubleshoot a few small errors that I came across in reconstructing the paper's tables.

#### Dayhoff et al. 1978

A big issue has always been poor quality for the scan of the paper that I found on the web.  The one I got through inter-library loan was no better (actually, it seems identical).  

There is a better [scan](http://compbio.berkeley.edu/class/c246/Reading/dayhoff-1978-apss.pdf) available now.

In particular, Fig 80 is very clear in the new version.  I put copies of these papers in the repo (``original papers`` directory), I  trust I won't get a takedown notice out of it.

#### original data

I found (many copies) of an old data file in my search of an archival hard drive.

* [dayhoff.corrected.changes.txt](data/dayhoff.corrected.changes.txt)

I don't recall the reason for the corrections, but the file consists of a lower triangular matrix of observed changes between amino acids in the proteins catalogued by Dayhoff and colleagues.

There is also appended what looks to be the outcome of a ``diff`` comparison of the corrected and original versions.
 
I decided to reformat this file with right-justified values, and spaces rather than tabs.  I wrote a short Python [script](script.py) for this.  The output was written to [dayhoff.fig80.txt](data/dayhoff.fig80.txt).

A last row was added listing the amino acids in order.

I noticed 2 fairly big errors in comparison with the cleaned up pdf of the paper that I found on the web.  I read the original as KR=417 but the clear copy is 477.  I also read the original as VP=150 but the clear copy shows 50.

Here is a [screenshot](screenshots/Dayhoff Fig 80.png) of Fig 80 from the new scan of the original paper.

### Redoing the analysis

#### Loading the data

The main Python script for working with the data is [dayhoff.py](dayhoff.py).  I'm going to work through the steps in the code and document them here.

The order of the amino acids may seem a little strange.  Here it is:

```
A R N D C Q E G H I L K M F P S T W Y V
```

They are in lexicographical order for the 3-letter abbreviations, but the symbols I'm using are 1-letter (so, for example, Arg=R comes at the second position).

The script saves that last row to [dayhoff.aaorder.txt](data/dayhoff.aaorder.txt), if the function ``save_order()`` is uncommented.

The first part of the code constructs a dictionary to hold the observed changes between amino acids in ``changeD``.  Changes are entered symmetrically, the value for **WD** is equal to the value for **DW**.

Next, the table is printed for proofreading.  Here is a [screenshot](screenshots/Screen Shot 2017-06-16 at 12.05.59 PM.png) of the output, which could be proofed against the input.

At this point we also tabulate the total number of changes observed involving the most observed changes (alanine, 3644 changes) compared to the least (tryptophan, 79 changes).  There is an important normalization that is needed before we can infer mutability, how changeable an amino acid is in evolutionary time.

#### Origin of Fig 80

What Dayhoff *et al.* did was to study 71 groups of closely related proteins known in 1978.  They observed a total of 1572 amino acid changes in these groups.  

They also inferred ancestral sequences:  i.e., given a tree:

```
(ACGH, DBGH), (ADIJ, CBIJ)
```

they infer ancestral sequences 

```
(ABGH, ABIJ)
```

The method for constructing the tree is not specified, but I believe they probably used parsimony.

They show this generic example:

```
ABGH -> ACGH and DBGH, infer A -> D and B -> C
ABIJ -> ADIJ and CBIJ, infer A -> C and B -> D
```

The first ancestor is not given but the changes are 

```
G <-> I and J <-> H.
```

From this analysis, they construct a table where AC and CA both have the value 1, etc.  Note that CD is 0, even though in this example the same position (2) has both C and D in different sequences, because it is inferred that this change was not direct, based on the tree.

The total of all the values in the dictionary is 31430 (next section).  Each change is counted twice, dividing by two we obtain 15715.  As described in the paper the total of all changes was 1572, the legend to Fig 80 explains that values were multiplied by 10.  The occurence of values in the table that are not multiples of 10 is explained by the comment

> In practice, some of the positions in the nodal sequences are blank [ambiguous]. For these, we have treated the changes statistically, distributing them among all observed alternatives.

#### Relative mutability

In order to say anything about the propensity of one amino acid to change, one into another (or simply to change at all), we need to know the frequencies of appearance of amino acids in the proteins Dayhoff et al studied.  This is provided in [Table 22](screenshots/Screen Shot 2017-06-16 at 12.26.44 PM.png).  I transcribed these as [dayhoff.frequencies.txt](data/dayhoff.frequencies.txt).  They were rewritten in row format for use with R ([dayhoff.frequencies.v2.txt](data/dayhoff.frequencies.v2.txt)).

Mutability is simply the ratio of observed changes to amino acid frequency, normalized to give ``Alanine = 100``.  Here is a [screenshot](screenshots/Screen Shot 2017-06-16 at 12.33.10 PM.png) of the output, which should be compared with [Table 21](screenshots/Dayhoff Table 21.png).

There are small discrepancies (compare 

| aa | theirs | ours |
| --- |:---:|:---:|
| N | 134 | 135 |
| S | 120 | 119 |
| E | 102 | 101 |
| T | 98 | 97 |
| H | 66 | 65 |
| W | 18 | 19 |

I'm not sure of the reason for this, but we've shown we can get very close to their values.  

From this point on, their values are used, in [dayhoff.mutabilities.txt](data/dayhoff.mutabilities.txt).

#### PAM1 matrix

The last step of this first half of the analysis is to generate the file [PAM1.txt](data/PAM1.txt), which should be compared with [Fig 82](screenshots/Dayhoff Fig 82.png) of the text.

> An element of this matrix, Mij, gives the probability that the amino acid in column j will be replaced by the amino acid in row i after a given evolutionary interval, in this case 1 PAM.
> 

As [described](screenshots/computation.png) in the paper, the off-diagonal entries in the PAM1 matrix (representing an amino acid not changing) are computed from the data in Fig 80.  

The original amino acids are specified in the columns (amino acid j), and the replacement amino acids are rows (amino acid i).

The first step in computation of the ij-th entry in PAM1 from the ij-th entry in Fig 80 is to form the ratio A<sub>ij</sub> / &Sigma;<sub>i</sub> A <sub>ij</sub>.  This is the ratio of observed changes involving amino acid i and j, divided by all changes involving amino acid j.

An additional factor is the mutability of amino acid j.  This seems kind of strange, since the mutability is &Sigma;<sub>i</sub> A <sub>ij</sub> divided by the frequency of j, f<sub>j</sub>.  In other words, &Sigma;<sub>i</sub> A <sub>ij</sub> cancels, leaving A<sub>ij</sub> / f<sub>j</sub>.

The last factor is a proportionality constant, &lambda;, which is the same for all columns.

> The quantity 100 X Zfi Mii gives the number of amino acids that will remain unchanged when a protein 100 links long, of average composition, is exposed to the evolutionary change represented by this matrix. This apparent evolutionary change depends upon the choice of X, in this case ,chosen so that this change is 1 mutation. Since there are almost no superimposed changes, this also represents 1 PAM of change. If h had been four times as large, the initial matrix would have represented 4 PAMs;  the discussion which follows would not be changed appreciably.

I didn't proof this carefully.  There are some differences which seem to be contained to the diagonal.  The on-diagonal values are supposed to be 10000 minus the sum of the rest of that column.

#### Continuing

You should now be able to follow the next [post](http://telliott99.blogspot.com/2008/08/pam-projecting-in-time.html) to generate the PAM250 matrix and so on.

If I get to it I will try to show how that part is done, but it is basically a matter of following the code in [Dayhoff.Rcode.txt](Dayhoff.Rcode.txt).