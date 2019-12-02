Here you can create a graph database with RNA transcripts, proteins, and genes from the Y chromosome of the dbSNP database.

The dbsnpY.py script downloads data from dbSNP, parses it, and writes it into CSV files.

The attached JSON files tell OrientDB how to read in each of the CSV files and what edges should be created.

The following questions arose and were answered during creation of this tutorial:

When making my schema, what should I start with at the center?
- start with the most specific, non-limiting thing, then branch out. For example, I started with RNA. While protein can be just as specific, ending at the end of the central dogma, sometimes there are ncRNAs that lead to no protein. Protein would therefore limit how many RNAs we are allowed to see. RNA:Protein has a 1:1 or 1:none relationship.

If I want a 1:many relationship, where does the foreign key belong?
- foreign key always belongs where there are many

When writing my CSV files, should I ensure that they are unique rows, or can I accomplish that in OrientDB?
- It is best to take care of this in the creation of CSVs.

Getting a java.lang.NumberFormatException? Do the first/id rows need to be numeric indices?
- Yes they must be ints

Why does SELECT * FROM Gene giving me only 20 nodes when there are 200 genes?
It is choosing the default/first cluster to show.

How do I change this so I can see all 200 genes?
- Go to the tool drop down menu and change limit

How is OrientDB automatically putting my data in clusters?
- https://orientdb.com/docs/2.2.x/Cluster-Selection.html

Can I assign them into clusters?
- https://orientdb.com/docs/2.2.x/Cluster-Selection.html#custom-cluster-selection-strategies

What is the max number of records in a cluster- why is it doing 20?
- https://orientdb.com/docs/2.2.x/Cluster-Selection.html
