# Nextflow Documentation Notes

[Nextflow](https://www.nextflow.io/)

[GitHub](https://github.com/nextflow-io/nextflow)

[Documentation](https://www.nextflow.io/docs/latest/index.html)

* Extends Unix pipes model.
* Dataflow programming model.
* Parallelisation implicitly defined by process input and output declarations.
* Continuous checkpointing of intermediate results => resume from last successful step.
* Supports Docker and Singularity containers.
* Ingegrates with GitHub.
* Executors for SGE, LSF, SLURM, PBS and HTCondor batch schedulers and for Kubernetes, Amazon AWS and Google Cloud platforms.

---

## Getting started

1. Check prerequisites

> Java 8 or later is required.

Install [OpenJDK](https://openjdk.java.net) 1.8:

* CentOS 7 users:

```console
$ sudo yum install -y openjdk-8-jdk-headless
```

* Ubuntu 18 users:

```console
$ sudo apt-get install -y openjdk-8-jdk-headless
```

Check install:

```console
$ javac -version
javac 1.8.0_242
$ java -version
openjdk version "1.8.0_242"
OpenJDK Runtime Environment (build 1.8.0_242-8u242-b08-0ubuntu3~18.04-b08)
OpenJDK 64-Bit Server VM (build 25.242-b08, mixed mode)
```

2. Set up

```console
$ curl -s https://get.nextflow.io | bash
```

This can also be installed via Bioconda.

3. Launch

```console
$ ./nextflow run hello
N E X T F L O W  ~  version 20.01.0
Pulling nextflow-io/hello ...
downloaded from https://github.com/nextflow-io/hello.git
Launching `nextflow-io/hello` [reverent_bose] - revision: 1d43afc0ec [master]
WARN: The use of `echo` method is deprecated
executor >  local (4)
[21/c2639b] process > sayHello [100%] 4 of 4 ?
Hello world!

Hola world!

Bonjour world!

Ciao world!
```

Note that content is pulled from https://github.com/nextflow-io/hello.

For future use:

```console
$ export PATH=$HOME/nextflow:$PATH
```

---

## Get started 

https://www.nextflow.io/docs/latest/getstarted.html

Write `tutorial.nf`:

```groovy
#!/usr/bin/env nextflow

params.str = 'Hello world!'

// Split a string into 6-character chunks.
// Write each one to a file with the prefix chunk_.
process splitLetters {
    output:
    file 'chunk_*' into letters

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

// Transform file contents to uppercase.
process convertToUpper {
    input:
    file x from letters.flatten()

    output:
    stdout result

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

// Emit results on "result channel" and print output using "view" operator.
result.view {it.trim()}
```

```console
$ nextflow run tutorial.nf 
N E X T F L O W  ~  version 20.01.0
Launching `tutorial.nf` [exotic_stonebraker] - revision: fc679f4860
executor >  local (3)
[fd/1c1343] process > splitLetters   [100%] 1 of 1 ?
[01/b40a76] process > convertToUpper [100%] 2 of 2 ?
WORLD!
HELLO
```

* `letters` is the common data between the processes.
* Process `convertToUpper` is executed in parallel, no guarantee of order.
* `fd/1c1343`: process ID and also prefix of process execution directory.

```console
$ ls -A work/fd/1c134388ee83af3f297ce9f562203f/
chunk_aa  .command.begin  .command.log  .command.run  .exitcode
chunk_ab  .command.err    .command.out  .command.sh
$ cat work/fd/1c134388ee83af3f297chunk_aa
Hello 
$ cat work/fd/1c134388ee83af3f297ce9f562hunk_ab
world!
```

Update `tutorial.nf`:

```groovy
// Reverse.
process convertToUpper {
    input:
    file x from letters

    output:
    stdout result

    """
    rev	$x
    """
}
```

Resume and rerun only processes that have changed, using cached results for unchanged processes:

```console
$ nextflow run tutorial.nf -resume
N E X T F L O W  ~  version 20.01.0
Launching `tutorial.nf` [backstabbing_shannon] - revision: a13f8c9089
executor >  local (1)
[fd/1c1343] process > splitLetters   [100%] 1 of 1, cached: 1 ?
[27/b9e6f7] process > convertToUpper [100%] 1 of 1 ?
olleH!dlrow
```

* `splitLetters` is skipped. process ID is unchanged, and results are retreived from cache.
* Default cache is `$PWD/work`

Define parameters using `params` prefix e.g. `params.str`.

Override parameters at command-line:

```console
$ nextflow run tutorial.nf --str 'Bonjour le monde'
N E X T F L O W  ~  version 20.01.0
Launching `tutorial.nf` [irreverent_stallman] - revision: a13f8c9089
executor >  local (2)
[43/503a6f] process > splitLetters   [100%] 1 of 1 ?
[ff/5f8883] process > convertToUpper [100%] 1 of 1 ?
uojnoBm el redno
```

---

## Basic concepts

https://www.nextflow.io/docs/latest/basic.html

Nextflow:

* Reactive workflow framework.
* Programming DSL.

Processes:

* Commands or scripts to be executed.
* Bash, Python, R etc.
* Isolated from each other.
* Executed independently.
* No common writable state.

Channels:

* Inputs and outputs.
* Inter-process communication via asynchronous FIFO queues.
* Implicitly defines pipeline execution flow.

```groovy
// Script parameters
params.query = "/some/data/sample.fa"
params.db = "/some/path/pdb"

db = file(params.db)
query_ch = Channel.fromPath(params.query)

process blastSearch {
    input:
    file query from query_ch

    output:
    file "top_hits.txt" into top_hits_ch

    """
    blastp -db $db -query $query -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits.txt
    """
}

process extractTopHits {
    input:
    file top_hits from top_hits_ch

    output:
    file "sequences.txt" into sequences_ch

    """
    blastdbcmd -db $db -entry_batch $top_hits > sequences.txt
    """
}
```

Channels and execution order:

* Execution order is determined by `blastSearch` and `extractTopHits` having a common `top_hits_ch` channel defined in their `output` and `input` respectively.
* These establish communication link between the process.
* `extractTopHits` waits for output and runs reactively.

Executor:

* Determines how scripts are run on a system.
* Local computer.
* HPC job submission: Open grid engine, Univa grid engine,Platform LSF, Linux SLURM, PBS Works, Torque, HTCondor
* Cloud: Amazon Web Services (AWS), Google Cloud Platform (GCP), Kubernetes

Scripting language:

* "designed to have a minimal learning curve, without having to pick up a new programming language"
* Domain Specific Language (DSL)
* Extension of Groovy programming language, a super-set of Java ("Python for Java")

Configuration:

* `nextflow.config`
* Pipeline configuration properties.
* Pipeline execution directory.
* Executors, environment variables, parameters etc.

```
process {
  executor = 'sge'
  queue = 'cn-el6'
}
```

---

## Nextflow scripting

https://www.nextflow.io/docs/latest/script.html

Execute any Groovy code or use any library for the JVM (Java) platform.

### Basic syntax

```groovy
println "Hello, World!"
```

Reuse variables:

```groovy
x = 1
x = new java.util.Date()
x = -3.1499392
x = false
x = "Hi"
```

Lists (`java.util.List`):

```groovy
myList = [1776, -1, 33, 99, 0, 928734928763]  // Zero-indexed.
println myList[0]
println myList.size()
(a, b, c) = [10, 20, 'foo']
println a
println b
println c
```

Maps (`java.util.Map`):

```groovy
scores = ["Brett":100, "Pete":"Did not finish", "Andrew":86.87934]
println scores["Pete"]
println scores.Pete
scores["Pete"] = 3
scores["Cedric"] = 120
```

Conditions:

```groovy
x = Math.random()
if (x < 0.5)
{
    println "You lost."
}
else
{
    println "You won!"
}
```

Strings:

```groovy
println "he said 'cheese' once"
println 'he said "cheese!" again'
a = "world"
print "hello " + a + "\n"

x = 'Hello'
println '$x + $y' // Literal string '$x + $y'.

foxtype = 'quick'
foxcolor = ['b', 'r', 'o', 'w', 'n']
println "The $foxtype ${foxcolor.join()} fox" // Variable substitution.

name = 'James'
text = '''
    hello there ${name}
    how are you today?
    '''

text = """
    hello there ${name}
    how are you today?
    """

myLongCmdline = """ blastp \
                -in $input_query \
                -out $output_file \
                -db $blast_database \
                -html
                """
result = myLongCmdline.execute().text
```

Closures:

```groovy
square = {it * it}
println square(3)            // 9
[1, 2, 3, 4].collect(square) // [1, 4, 9, 16]

printMapClosure = {key, value ->
    println "$key = $value"
}
["Yue": "Wu", "Mark": "Williams", "Sudha": "Kumari"].each(printMapClosure)

myMap = ["China": 1 , "India": 2, "USA": 3]
result = 0
myMap.keySet().each({result += myMap[it]})
println result
```

* `it` is implicit variable for value passed.
* `collect` applies closure to each item, creating a new list.
* `each` applies two-argument closure to each key-value, creating a new map.
* Similar to Python list/dictionary comprehension.
* Can access variables within the scope in which they are defined.

### Regular expressions

Regular expressions:

```groovy
assert 'foo' =~ /foo/    // TRUE
assert 'foobar' =~ /foo/ // TRUE

assert 'foo' ==~ /foo/    // TRUE
assert 'foobar' ==~ /foo/ // FALSE

x = ~/abc/
println x.class            // java.util.regex.Pattern
y = 'some string' =~ /abc/
println y.class            // java.util.regex.Matcher

programVersion = '2.7.3-beta'
m = programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/
assert m[0] ==  ['2.7.3-beta', '2', '7', '3', 'beta'] // 1st regexp match
assert m[0][1] == '2'    // 1st group in regexp
assert m[0][2] == '7'    // 2nd group in regexp
assert m[0][3] == '3'
assert m[0][4] == 'beta'

(full, major, minor, patch, flavor) = (programVersion =~ /(\d+)\.(\d+)\.(\d+)-?(.+)/)[0]
println full   // 2.7.3-beta
println major  // 2
println minor  // 7
println patch  // 3
println flavor // beta
```

[java.util.regex.Pattern](https://docs.oracle.com/javase/7/docs/api/java/util/regex/Pattern.html):

* `~/pattern/`
* `=~`: pattern occurs anywhere in a string?
* `==~`: string matches pattern exactly?
* `(...)`: capture group.
* `\d+`: 1+ digits.
* `\.`: dot character.
* `-?`: zero or one `-`.
* `.+`: 1+ characters.

String replacement and substring removal:

```groovy
x = "colour".replaceFirst(/ou/, "o")            // 'color'
x = "cheesecheese".replaceAll(/cheese/, "nice") // 'nicenice'

wordStartsWithGr = ~/(?i)\s+Gr\w+/
x = 'Hello Groovy world!' - wordStartsWithGr // 'Hello world!'
x = 'Hi Grails users' - wordStartsWithGr     // 'Hi users'
```

* `(?i)`: case-insensitive matching.
* `\s+`: 1+ more whitespace characters.
* `\w+`: 1+ more word characters.

```groovy
x = 'Remove first matching 5 letter word' - ~/\b\w{5}\b/) // 'Remove  matching 5 letter word'
x = 'Line contains 20 characters' - ~/\d+\s+/ // 'Line contains characters'
```

* `\b`: word boundary.
* `\w{5}`: 5 word characters.
* `\d+`: one or more digits.
* `\s+`: one or more whitespace characters.

### Files (`java.nio.file.Path`) and I/O

```groovy
myFile = file('some/path/to/my_file.file')
listOfFiles = file('some/path/*.fa')
listWithHidden = file('some/path/*.fa', hidden: true)
```

glob patterns:

* `*`: 0+ characters.
* `**`: 0+ characters, and cross directory boundaries.
* `?`: single character.
* `[abc]`: match any single character.
* `{sun,moon,stars}`: match any word.
* If glob pattern provided, `file` returns list of matching files or `[]`.

Read/write all content to/from a single variable:

```groovy
print myFile.text             // Load and print.
myFile.text = 'Hello world!'  // Save, overwrite existing content.
myFile.append('Add a line\n')
myFile << 'Add another line\n'

binaryContent = myFile.bytes  // Byte array
myFile.bytes = binaryBuffer
```

Read all content to a list:

```groovy
myFile = file('some/my_file.txt')
allLines  = myFile.readLines()
for(line: allLines) {
    println line
}

file('some/my_file.txt')
    .readLines()
    .each {println it}
```

Read content line-by-line:

```groovy
count = 0
myFile.eachLine { str ->
        println "line ${count++}: $str"
}
```

Read files:

```groovy
myReader = myFile.newReader()
String line
while(line = myReader.readLine()) {
    println line
}
myReader.close()

myFile.withReader {
    String line
    while(line = myReader.readLine()) {
        println line
    }
}
```

* Text:
  - `newReader`, `withReader`, `java.io.Reader`
  - Read single characters, lines or arrays of characters.
* Binary: `newInputStream`, `withInputStream`, `java.io.InputStream`

Write files:

```groovy
sourceFile.withReader { source ->
    targetFile.withWriter { target ->
        String line
        while(line=source.readLine()) {
            target << line.replaceAll('U','X')
        }
    }
}
```

* Text: `newWriter`, `withWriter`, `java.io.Writer`
* Formatted text:`newPrintWriter`, `withPrintWriter`, `java.io.PrintWriter`
* Binary: `newOutputStream`, `withOutputStream`,  `java.io.OutputStream`

List directory:

```groovy
myDir = file('any/path')
allFiles = myDir.list()
for(def file : allFiles) {
    println file
}

myDir.eachFile { item ->
    if (item.isFile())
    {
        println "${item.getName()} - size: ${item.size()}"
    }
    else if (item.isDirectory())
    {
        println "${item.getName()} - DIR"
    }
}
```

`list`, `listFiles`, `eachDir`, `eachFileMatch`, `eachDirMatch`, `eachFileRecurse`, `eachDirRecurse`

Note: for `def` in `for(def file : allFiles)`: Stackoverflow's [What is the difference between defining variables using def and without?](https://stackoverflow.com/questions/39514795/what-is-the-difference-between-defining-variables-using-def-and-without/39546627) comments that:

> When you assign a value to a variable without a "def" or other type ... it's added to the "binding", the global variables for the script. That means it can be accessed from all functions within the script. It's a lot like if you had the variable defined at the top of the script.

Directory and file manipulation:

```groovy
result = myDir.mkdir()  // TRUE or FALSE
result = myDir.mkdirs() // Create non-existent parent directories.

myFile.copyTo('new_name.txt')
myDir.copyTo('/some/new/path')

myFile.moveTo('/another/path/new_file.txt')
myDir = file('/any/dir_a')
myDir.moveTo('/any/dir_b') // If dir_b exists then dir_a is moved into it
                           // Otherwise, dir_a is renamed to dir_b

myFile.renameTo('new_file_name.txt')

result = myFile.delete() // TRUE or FALSE

myDir.delete()    // Delete empty directory only
myDir.deleteDir()
```

* `getName`: `/some/path/file.txt` => `file.txt`
* `getBaseName`: `/some/path/file.tar.gz` => `file.tar`
* `getSimpleName`: `/some/path/file.tar.gz` => `file`
* `getExtension`: `/some/path/file.txt` => `txt`
* `getParent`: `/some/path/file.txt` => `/some/path`
* `exists`, `isEmpty`, `isFile`, `isDirectory`

### HTTP/FTP files

```groovy
pdb = file('http://files.rcsb.org/header/5FID.pdb')
```

### Counting records

```groovy
fa = file('/home/ubuntu/riboviz/vignette/input/yeast_rRNA_R64-1-1.fa')
println "${fa}"
println "Lines: ${fa.countLines()}"
println "FASTA records: ${fa.countFasta()}"
fq = file('/home/ubuntu/riboviz/data/simdata/umi5_umi3.fastq')
println "${fa}"
println "Lines: ${fq.countLines()}"
println "FASTQ records: ${fq.countFastq()}"
fq = file('/home/ubuntu/riboviz/vignette/input/SRR1042855_s1mi.fastq.gz')
println "${fq}"
println "Lines: ${fq.countLines()}"
println "FASTQ records: ${fq.countFastq()}"
```
```console
/home/ubuntu/riboviz/vignette/input/yeast_rRNA_R64-1-1.fa 
Lines: 205
FASTA records: 12
/home/ubuntu/riboviz/data/simdata/umi5_umi3.fastq 
Lines: 36
FASTQ records: 9
/home/ubuntu/riboviz/vignette/input/SRR1042855_s1mi.fastq.gz 
Lines: 3854284
FASTQ records: 963571
```

---

## Processes

https://www.nextflow.io/docs/latest/process.html

```groovy
process < name > {
   [ directives ]
   input:
    < process inputs >
   output:
    < process outputs >
   when:
    < condition >
   [script|shell|exec]:
   < user script to be executed >
}
```
```groovy
process sayHello {
    """
    echo 'Hello world!' > file
    """
}
```

### Script

Commands to be executed via bash.

```groovy
process doMoreThings {
    """
    blastp -db $db -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    blastdbcmd -db $db -entry_batch top_hits > sequences
    """
}
```

* `'...'` strings are literal.
* `"..."` strings support substitutions.
* Take care if variable substitutions are to be evaluated by Nextflow or by bash.

Environment variables:

```groovy
process printPath {
    '''
    echo The path is: $PATH
    '''
}
```

* Cannot use Nextflow variables.
* Escape environment variables.

```
process doOtherThings {
    """
    blastp -db \$DB -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n $MAX | cut -f 2 > top_hits
    blastdbcmd -db \$DB -entry_batch top_hits > sequences
    """
}
```

Use `#!/usr/bin/env <INTERPRETER>` for cross-platform portability.

Python:

```groovy
process pythonExample {
    output:
    stdout python_result

    """
    #!/usr/bin/env python

    x = 'Hello'
    y = 'Python!'
    print("%s - %s" % (x,y))
    """
python_result.view { it.trim() }
```

```groovy
x = 'Hello'
y = 'Python!'
process pythonExample {
    output:
    stdout python_result

    """
    #!/usr/bin/env python

    print('${x}')
    print("%s - %s" % ("${x}", "${y}"))
    """
}
python_result.view { it.trim() }
```

R:

```groovy
process rExample {
    output:
    stdout r_result

    """
    #!/usr/bin/env Rscript

    x <- 'Hello'
    y <- 'R!'
    print(paste(x, y))
    """
}
r_result.view { it.trim() }
```
```groovy
x = 'Hello'
z = 'R!'
process rExample {
    output:
    stdout r_result

    """
    #!/usr/bin/env Rscript

    print(paste("${x}", "${z}"))
    """
}
r_result.view { it.trim() }
```

`script` blocks and conditional execution:

```groovy
seq_to_align = ...
mode = 'tcoffee'

process align {
    input:
    file seq_to_aln from sequences

    script:
    if (mode == 'tcoffee')
        """
        t_coffee -in $seq_to_aln > out_file
        """
    else if (mode == 'mafft')
        """
        mafft --anysymbol --parttree --quiet $seq_to_aln > out_file
        """
    else if (mode == 'clustalo')
        """
        clustalo -i $seq_to_aln -o out_file
        """
    else
        error "Invalid alignment mode: ${mode}"
}
```

* `script` block must evaluate to the script to be executed.
* Supports `if`, `switch` etc.

Templates:

```groovy
process template_example {
    input:
    val STR from 'this', 'that'

    script:
    template 'my_script.sh'
}
```

* `templates/my_script.sh`:

```console
#!/bin/bash
echo "process started at `date`"
echo $STR
:
echo "process completed"
```

* If located elsewhere then specify absolute path in `template`.
* `$STR` evaluation:
  - Bash variable, if `my_script.sh` is run standalone.
  - Nextflow variable placeholder, if run as a template.

`shell` blocks:

```groovy
process myTask {
    input:
    val str from 'Hello', 'Hola', 'Bonjour'

    shell:
    '''
    echo User $USER says !{str}
    '''
}

process doOtherThings {
    '''
    blastp -db $DB -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n !{MAX} | cut -f 2 > top_hits
    blastdbcmd -db $DB -entry_batch top_hits > sequences
    '''
}
```

* Use single-quoted strings only.
* Use `!{...}` for Nextflow variable substitions.
* Can use with `template`.

`exec` blocks and native Groovy execution:

```groovy
x = Channel.from('a', 'b', 'c')

process simpleSum {
    input:
    val x

    exec:
    println "Hello Mr. $x"
}
```

### Inputs

```groovy
input:
    <input qualifier> <input name> [from <source channel>] [attributes]
```

* If `input name` is same as channel name then `from` can be omitted.

`val`:

```groovy
num = Channel.from(1, 2, 3)
process basicExample {
    input:
    val x from num

    "echo process job $x"
}

process alternateBasicExample {
    input:
    val num

    "echo process job $num"
}
```

`file`:

```groovy
proteins = Channel.fromPath('/some/path/*.fa')

process blastThemAll {
    input:
    file query_file from proteins

    "blastp -query ${query_file} -db nr"
}

process alternateBlastThemAll {
    input:
    file proteins

    "blastp -query $proteins -db nr"
}
```

Stage different files with same file name:

```groovy
input:
    file query_file name 'query.fa' from proteins
```

* or:

```groovy
input:
    file 'query.fa' from proteins
```

```groovy
proteins = Channel.fromPath('/some/path/*.fa')
process blastThemAll {
    input:
    file 'query.fa' from proteins

    "blastp -query query.fa -db nr"
}
```

* Each file is staged with name `query.fa` in a different execution context.

Multiple input files and controlling file naming:

```groovy
fasta = Channel.fromPath("/home/ubuntu/riboviz/vignette/input/*.fa").buffer(size:2)

process blastThemAll {
    input:
    file 'seq' from fasta

    output:
    stdout blast_result

    """
    echo seq*
    wc -l seq*
    """
}
```
```console
seq1 seq2
   205 seq1
  1851 seq2
  2056 total
```

* Determines how file is named in work directory.
* Target input file name can include globs.

```groovy
fasta = Channel.fromPath("/some/path/*.fa").buffer(size:3)

process blastThemAll {
    input:
    file 'seq?.fa' from fasta

    "cat seq1.fa seq2.fa seq3.fa"
}
```

* To preserve original file names use patterns or variable identifiers, not '<FILENAME>'.

Dynamic input file names:

* Use variable name in file name:

```groovy
process simpleCount {
    input:
    val x from species
    file "${x}.fa" from genomes

    """
    cat ${x}.fa | grep '>'
    """
}
```

`path`:

* Backwards compatible with `file`.
* Interpret strings as paths,
* Must be prefixed with `/` or a URI protocol.

```groovy
process foo {
    input:
    path x from '/some/data/file.txt'

    """
    your_command --in $x
    """
}
```

* Control file naming:

```groovy
process foo {
    input:
    path x, stageAs: 'data.txt' from '/some/data/file.txt'

    """
    your_command --in data.txt
    """
}
```

`stdin`:

```groovy
str = Channel.from('hello', 'hola', 'bonjour', 'ciao').map { it+'\n' }
process printAll {
    input:
    stdin str

    """
    cat -
    """
}
```

`env`:

```groovy
str = Channel.from('hello', 'hola', 'bonjour', 'ciao')
process printEnv {
    input:
    env HELLO from str

    '''
    echo $HELLO world!
    '''
}
```

`tuple`:

```groovy
values = Channel.of([1, 'alpha'], [2, 'beta'], [3, 'delta'])
process tupleExample {
    input:
    tuple val(x), val(y) from values

    output:
    stdout tuple_result

    """
    echo Processing $x $y
    """
}
```

* In the following, `cat` receives its input from `file`.

```groovy
process tupleExample {
    input:
    tuple val(x), file('latin.txt') from values

    output:
    stdout tuple_result

    """
    echo Processing $x
    cat - latin.txt > copy
    """
}

process alternateTupleExample {
    input:
    tuple x, 'latin.txt' from values

    """
    echo Processing $x
    cat - latin.txt > copy
    """
}
```

Permutations, `each` and input repeaters:

```groovy
sequences = Channel.fromPath('*.fa')
methods = ['regular', 'expresso', 'psicoffee']
process alignSequences {
    input:
    file seq from sequences
    each mode from methods

    """
    t_coffee -in $seq -mode $mode > result
    """
}
```

```groovy
sequences = Channel.fromPath('*.fa')
methods = ['regular', 'expresso']
libraries = [ file('PQ001.lib'), file('PQ002.lib'), file('PQ003.lib') ]
process alignSequences {
    input:
    file seq from sequences
    each mode from methods
    each file(lib) from libraries

    """
    t_coffee -in $seq -mode $mode -lib $lib > result
    """
}
```

Multiple input channels:

* Process waits for input on all input channels.
* Consumes values serially.
* Execution stops if any channel has no more values.

```groovy
process foo {
    echo true

    input:
    val x from Channel.from(1,2)
    val y from Channel.from('a','b','c')

    script:
    """
    echo $x and $y
    """
}
```

* Executes using `1` and `a` then `2` abd `b`. `c` is discarded.
* Value/singleton channel can be read indefinitely:

```groovy
process bar {
    echo true

    input:
    val x from Channel.value(1)
    val y from Channel.from('a','b','c')
  
    script:
    """
    echo $x and $y
    """
}
```

### Outputs

```groovy
output:
    <output qualifier> <output name> [into <target channel>[,channel,..]] [attribute [,..]]
```

* If output name is same as channel name then `into` can be omitted.

`val`:

```groovy
methods = ['prot','dna', 'rna']
process foo {
    input:
    val x from methods

    output:
    val x into receiver
    """
    echo $x > file
    """
}
receiver.println { "Received: $it" }
```

```groovy
process foo {
    input:
    file fasta from 'dummy'

    output:
    val x into var_channel
    val 'BB11' into str_channel
    val "${fasta.baseName}.out" into exp_channel

    script:
    x = fasta.name
    """
    cat $x > file
    """
}
```

`file`:

```groovy
process randomNum {
    output:
    file 'result.txt' into numbers

    '''
    echo $RANDOM > result.txt
    '''
}
numbers.subscribe { println "Received: " + it.text }
```

Multiple output files:

* Capture multiple output files into a list and output them.
* Output file names with `?` or `*`.
* Matching input files are not matched to the output but they are copied to downstream tasks.

```groovy
process splitLetters {
    output:
    file 'chunk_*' into letters

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
letters
    .flatMap()
    .subscribe { println "File: ${it.name} => ${it.text}" }
```
```console
File: chunk_aa => H
File: chunk_ab => o
File: chunk_ac => l
File: chunk_ad => a
```

* `printf "Hola" | split -b 1 - chunk_` creates `chunk_aa`, `chunk_ab`, `chunk_ac`, `chunk_ad` with `H`, `o`, `l`, `a`.
* `mode flatten` (deprecated as of 19.10.0):

```groovy
process splitLetters {
    output:
    file 'chunk_*' into letters mode flatten

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
letters.subscribe { println "File: ${it.name} => ${it.text}" }
```

Dynamic output file names:

```groovy
process align {
    input:
    val x from species
    file seq from sequences

    output:
    file "${x}.aln" into genomes

    """
    t_coffee -in $seq > ${x}.aln
    """
}
```

> in most cases, you don't need to take care of naming output files, because each task is executed in its own unique temporary directory, so files produced by different tasks can never override each other.
>
> meta-data can be associated with outputs by using the tuple output qualifier, instead of including them in the output file name.
>
> use of output files with static names over dynamic ones is preferable whenever possible, because it will result in a simpler and more portable code

`path`:

* Backwards compatible with `file`.
* `file` interprets `:` as path separator.
* `path` interprets `:` as file name character.

`stdout`:

```groovy
process echoSomething {
    output:
    stdout channel

    """
    echo Hello world!
    """
}
channel.subscribe { print "I say..  $it" }
```

`env`:

```groovy
process myTask {
    output:
    env FOO into target

    script:
    '''
    FOO=$(ls -la)
    '''
}
target.view { "directory content: $it" }
```

`tuple`:

```groovy
query_ch = Channel.fromPath '*.fa'
species_ch = Channel.from 'human', 'cow', 'horse'
process blast {
    input:
    val species from query_ch
    file query from species_ch

    output:
    tuple val(species), file('result') into blastOuts

    script:
    """
    blast -db nr -query $query > result
    """
}
```

* Short form:

``` groovy
    output:
    tuple species, 'result' into blastOuts
```

Optional output:

```groovy
    output:
    file("output.txt") optional true into outChannel
```

### When and conditionals

```groovy
process find {
    input:
    file proteins
    val type from dbtype

    when:
    proteins.name =~ /^BB11.*/ && type == 'nr'

    script:
    """
    blastp -query $proteins -db nr
    """
}
```

### Directives

Optional settings that affect process execution.

```groovy
    name value [, value2 [,..]]
```

Myriad directives, examples include...

`beforeScript` / `afterScript`:

```groovy
process foo {
    beforeScript 'source /cluster/bin/setup'

    """
    echo bar
    """
}
```

`cache true`:

* Store process results in local cache.
* Enabled by default.
* Run Nextflow with `resume` and same inputs to produce stored data as actual results.
* Generates unique key based on process script and inputs.
* Can be disabled on a process-specific basis.

`conda 'bwa=0.7.15 ...'`: define process dependencies using Conda.

`container 'dockerbox:tag'`: execute process script in a Docker containe.

`cpus <N>`: number of CPUs required by process task, so these can be reserved.

`echo true`: forward all command stdout to standard output.

`errorStrategy: terminate|finish|ignore|retry`

`publishDir '/.../'`:

* Publish process output files to specific folder.

```groovy
process foo {
    publishDir '/data/chunks'

    output:
    file 'chunk_*' into letters

    '''
    printf 'Hola' | split -b 1 - chunk_
    '''
}
```

`label`: annotate process with mnemonic identifier.

`storeDir '/db/genomes'`: define directory as permanent cache for process results.

`stageInMode copy|link|symlink|rellink`: input file staging.

`stageOutMode copy|move|rsync`: output file staging.

`tag` associate process execution with custom label.

`validExitStatus 0,1,2`: valid exit status for process.

Dynamic directives:

```groovy
process foo {
    executor 'sge'
    queue { entries > 100 ? 'long' : 'short' }

    input:
    set entries, file(x) from data

    script:
    """
    < your job here >
    """
}
```

* Access dynamic directive value using `${task.<DIRECTIVE>}`.

Dynamic computing resources and `task.attempt`:

```groovy
process foo {
    memory { 2.GB * task.attempt }
    time { 1.hour * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    script:
    <your job here>
}
```

---

## Channels

https://www.nextflow.io/docs/latest/channel.html

Behaviour:

* Sending data is asynchronous, completes immediately.
* Receiving data is blocking, stopsreceiving process until data has arrived.

### Channel types

Queue channel:

* Non-blocking unidirectional FIFO queue.
* Factory method e.g. `from`, `fromPath` etc.
* Channel operator e.g. `map`, flatMap`.
* `output` declaration `into` clause e.g. `file 'chunk_*' into letters`. 
* Same queue channel cannot be used more than one time as process output and more than one time as process input - use `into` operator to create copies.

Value channel:

* AKA singleton channel.
* Read unlimited times without consuming its content.
* Factory method.
* Operator returning single value e.g. `first`, `last`, `collect`, `count`, `min`, `max`, `reduce`, `sum`.
* Implicitly created if single value specified in `from` clause.
* Implicitly created for output if process has only value channels as inputs.

```groovy
process foo {
    input:
    val x from 1

    output:
    file 'x.txt' into result

    """
    echo $x > x.txt
    """
}
result.view { it }
```
```console
/home/ubuntu/workflows/nextflow/work/29/4cc7bc1607fa908ef027a59f7022e0/x.txt
$ cat /home/ubuntu/workflows/nextflow/work/29/4cc7bc1607fa908ef027a59f7022e0/x.txt
1
```

### Channel factory

`of`:

```groovy
Channel.of(1, 3, 5, 7) // 1 3 5 7
Channel.of(1..23, 'X', 'Y') // 1 2 ... 23 X Y
```

`from`:

* deprecated, use `of` or `fromList`.
* Identical:

```groovy
Channel.from(1, 3, 5, 7, 9)
Channel.from([1, 3, 5, 7, 9])
```

* Not identical:

```groovy
Channel.from(1, 3, 5, 7, 9)
Channel.from([1, 3], [5, 7], [9]])
```

`value`:

```
expl1 = Channel.value()                // Empty variable
expl2 = Channel.value('Hello there')
expl3 = Channel.value([1, 2, 3, 4, 5]) // List value
```

`fromList`:

```groovy
Channel.fromList(['a', 'b', 'c', 'd']) // a b c d
```

`fromPath`:

```groovy
myFileChannel = Channel.fromPath('/data/some/bigfile.txt')

* Emits one path per match for globs.

```groovy
myFileChannel = Channel.fromPath('/data/big/*.txt')
files = Channel.fromPath('data/**.fa')
moreFiles = Channel.fromPath('data/**/*.fa')
pairFiles = Channel.fromPath('data/file_{1,2}.fq')
```

* Hidden files:

```groovy
expl1 = Channel.fromPath('/path/.*')
expl2 = Channel.fromPath('/path/.*.fa')
expl3 = Channel.fromPath('/path/*', hidden: true)
```

* File types:

```groovy
myFileChannel = Channel.fromPath('/path/*b', type: 'dir')
myFileChannel = Channel.fromPath('/path/a*', type: 'any')
```

* Lists of paths:

```groovy
Channel.fromPath(['/some/path/*.fq', '/other/path/*.fastq'])
```

`fromFilePairs`:

* Pairs of grouping key of match and file list.

```groovy
Channel.fromFilePairs('/my/data/SRR*_{1,2}.fastq')
```
```console
[SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
[SRR493367, [/my/data/SRR493367_1.fastq, /my/data/SRR493367_2.fastq]]
[SRR493368, [/my/data/SRR493368_1.fastq, /my/data/SRR493368_2.fastq]]
[SRR493369, [/my/data/SRR493369_1.fastq, /my/data/SRR493369_2.fastq]]
[SRR493370, [/my/data/SRR493370_1.fastq, /my/data/SRR493370_2.fastq]]
[SRR493371, [/my/data/SRR493371_1.fastq, /my/data/SRR493371_2.fastq]]
```

* Custom grouping key e.g. extensions and files:

```groovy
Channel.fromFilePairs('/some/data/*', size: -1) { file -> file.extension }
       .println { ext, files -> "Files with the extension $ext are $files" }
```

* Lists of paths:

```groovy
Channel.fromFilePairs(['/some/data/SRR*_{1,2}.fastq', '/other/data/QFF*_{1,2}.fastq'])
```

`fromSRA`:

* Query NCBI SRA database and emit FASTQ files matching specified criteria.
* Supports any query supported by NCBI ESearch API.
* Specify API key:
  - `api_key` optional parameter.
  - `NCBI_API_KEY` environment variable.

```groovy
Channel.fromSRA('SRP043510')
```
```console
[SRR1448794, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/004/SRR1448794/SRR1448794.fastq.gz]
[SRR1448795, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/005/SRR1448795/SRR1448795.fastq.gz]
[SRR1448792, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR144/002/SRR1448792/SRR1448792.fastq.gz]
...

ids = ['ERR908507', 'ERR908506', 'ERR908505']
Channel.fromSRA(ids)
[ERR908507, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908507/ERR908507_2.fastq.gz]]
[ERR908506, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908506/ERR908506_2.fastq.gz]]
[ERR908505, [ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_1.fastq.gz, ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR908/ERR908505/ERR908505_2.fastq.gz]]
```

`watchPath`:

* Emit file meeting specific condition (`create`, `modify`, `delete`).
* Yields non-terminating pipeline.

```groovy
Channel.watchPath('/path/*.fa')
       .subscribe { println "Fasta file: $it" }
Channel.watchPath('/path/*.fa', 'create,modify')
       .subscribe { println "File created or modified: $it" }
```

`empty`:

* Does not emit any value.

### Binding values

```groovy
myChannel = Channel.create()
myChannel.bind('Hello world')
```
```groovy
myChannel = Channel.create()
myChannel << 'Hello world'
```

### Observing events

`subscribe`:

* Execute user-defined function (closure) every time a value is emitted.

```groovy
source = Channel.from('alpha', 'beta', 'delta')
source.subscribe { println "Got: $it" }
```
```console
Got: alpha
Got: beta
Got: delta
```

```groovy
Channel.from('alpha', 'beta', 'lambda')
       .subscribe { String str ->
           println "Got: ${str}; len: ${str.size()}"
        }
```
```console
Got: alpha; len: 5
Got: beta; len: 4
Got: lambda; len: 6
```

* `onNext`: whenever value is emitted, equivalent to example above.
* `onComplete`: after last value is emitted.
* `onError`: when exception is raised in `onNext`.

```groovy
Channel.from(1, 2, 3)
       .subscribe onNext: { println it }, onComplete: { println 'Done' }
```
```console
1
2
3
Done
```

---

## Operators

https://www.nextflow.io/docs/latest/operator.html

Connect channels and transform values.

### Filtering operators

```groovy
Channel.from('a', 'b', 'aa', 'bc', 3, 4.5).filter(~/^a.*/) // a aa
Channel.from('a', 'b', 'aa', 'bc', 3, 4.5).filter(Number) // 3 4.5
Channel.from(1, 2, 3, 4, 5).filter {it % 2 == 1} // equivalent to filter({it.toString().size() == 1}). 1 3 5
```
```groovy
Channel.from(1, 1, 1, 5, 7, 7, 7, 3, 3).unique() // 1 5 7 3
Channel.from(1, 3, 4, 5).unique { it % 2 } // 1 4
```
```groovy
Channel.from(1, 1, 2, 2, 2, 3, 1, 1, 2, 2, 3).distinct() // 1 2 3 1 2 3
Channel.from(1, 1, 2, 2, 2, 3, 1, 1, 2, 4, 6).distinct { it % 2 } // 1 2 3 2
          // 1  1  0  0  0  1  1  1  0  0  0
    
```
```groovy
Channel.from(1, 2, 3).first() // 1
Channel.from(1, 2, 'a', 'b', 3).first(String) // a
Channel.from('a', 'aa', 'aaa').first(~/aa.*/) // aa
Channel.from(1, 2, 3, 4, 5).first { it > 3 } // 4
```
```groovy
Channel.from(1..100 ).randomSample(10) // 10 random numbers between 1..100
Channel.from(1..100).randomSample(10, 234) // Seed with 234.
```
```groovy
Channel.from(1, 2, 3, 4, 5, 6).take(3) // 1 2 3
```

* Use -1 for all values.

```groovy
Channel.from(1, 2, 3, 4, 5, 6).last() // 6
```
```groovy
Channel.from(3, 2, 1, 5, 1, 5).until{ it==5 } // 3 2 1
```

### Transforming operators

```groovy
Channel.from(1, 2, 3, 4, 5).map {it * it} // 1 4 9 16 25
```
```groovy
numbers = Channel.from(1, 2, 3)
results = numbers.flatMap { n -> [ n*2, n*3 ] } // Flattens any lists, 2 3 4 6 6 9
Channel.from(1, 2, 3).flatMap { it -> [ number: it, square: it*it ] } // number: 1 square: 1 ....
```

`reduce` takes result of previous `reduce` (accumulator) and emitted item:

```groovy
Channel.from(1, 2, 3, 4, 5)
       .reduce { a, b -> println "a: $a b: $b"; return a+b }
       .subscribe { println "result = $it" }
```
```console
a: 1    b: 2
a: 3    b: 3
a: 6    b: 4
a: 10   b: 5
result = 15
```
```groovy
myChannel.reduce(seedValue) {  a, b -> ... }
```
```groovy
Channel.from('hello', 'ciao', 'hola', 'hi', 'bonjour')
       .groupBy { String str -> str[0] }
       .subscribe { println it }
```
```console
[ b:['bonjour'], c:['ciao'], h:['hello','hola','hi'] ]
```

`groupTuple` 

* Collect tuples or lists and group those with common keys.
* Transform list[(K, V, W, ..)] into list(K, list(V), list(W), ..)

```groovy
Channel.from([1, 'A'], [1, 'B'], [2, 'C'], [3, 'B'], [1, 'C'], [2, 'A'], [3, 'D'])
       .groupTuple()
       .subscribe { println it }
```
```console
[1, [A, B, C]]
[2, [C, A]]
[3, [B, D]]
```

`buffer`:

* `buffer(openingCondition, closingCondition)`:

```groovy
Channel.from(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2).buffer(2, 4) // [2, 3, 4] [2, 3, 4]
```

* `buffer(size: n`)`:

 transform the source channel in such a way that it emits tuples made up of n elements. An incomplete tuple is discarded. For example:

```groovy
Channel.from(1, 2, 3, 1, 2, 3, 1).buffer(size: 2) // [1, 2] [3, 1] [2, 3]
Channel.from(1, 2, 3, 1, 2, 3, 1).buffer(size: 2, remainder: true) // [1, 2] [3, 1] [2, 3] [1]
```

* `buffer(size: n, skip: m)`:

```groovy
Channel.from(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2).buffer(size:3, skip:2) // [3, 4, 5] [3, 4, 5]
```

`collate`:

* Collate into tuples of N items.

```groovy
Channel.from(1, 2, 3, 1, 2, 3, 1).collate(3) // [1, 2, 3] [1, 2, 3] [1]
Channel.from(1, 2, 3, 1, 2, 3, 1).collate(3, false) // [1, 2, 3] [1, 2, 3]
Channel.from(1, 2, 3, 4).collate(3, 1) // [1, 2, 3] [2, 3, 4] [3, 4] [4]
```

`collect`:

```groovy
Channel.from(1, 2, 3, 4).collect() // [1, 2, 3, 4]
Channel.from('hello', 'ciao', 'bonjour').collect { it.length() } // [5, 4, 7]
```

`flatten`:

```groovy
Channel.from([1, [2, 3]], 4, [5, [6]]).flatten() // 1 2 3 4 5 6
```

`toList`:

```groovy
Channel.from(1, 2, 3, 4).toList() // [1,2,3,4]
```

`toSortedList`:

```groovy
Channel.from(3, 2, 1, 4).toSortedList() // [1,2,3,4]
Channel
    .from(["homer", 5], ["bart", 2], ["lisa", 10], ["marge", 3], ["maggie", 7])
    .toSortedList({ a, b -> b[1] <=> a[1] }) // [[lisa, 10], [maggie, 7], [homer, 5], [marge, 3], [bart, 2]]
```

* `A <=> B` return -1 if A < B, 0 if A == B, 1 if A > B.

`transpose`:

```groovy
Channel.from([
   ['a', ['p', 'q'], ['u','v'] ],
   ['b', ['s', 't'], ['x','y'] ]
   ])
   .transpose()
  .println()
```
```console
[a, p, u]
[a, q, v]
[b, s, x]
[b, t, y]
```

### Splitting operators

`splitCsv`:

```groovy
Channel.from('alpha,beta,gamma\n10,20,30\n70,80,90')
       .splitCsv()
       .subscribe { row ->
          println "${row[0]} - ${row[1]} - ${row[2]}"
       }
```
```console
alpha - beta - gamma
10 - 20 - 30
70 - 80 - 90
```
```groovy
Channel.from('alpha,beta,gamma\n10,20,30\n70,80,90')
       .splitCsv(header: true)
       .subscribe { row ->
          println "${row.alpha} - ${row.beta} - ${row.gamma}"
       }
```
```console
10 - 20 - 30
70 - 80 - 90
```
```groovy
Channel.from('alpha,beta,gamma\n10,20,30\n70,80,90')
       .splitCsv(header: ['col1', 'col2', 'col3'], skip: 1 )
       .subscribe { row ->
          println "${row.col1} - ${row.col2} - ${row.col3}"
       }
```

`splitFasta`:

```groovy
Channel.fromPath('misc/sample.fa')
       .splitFasta(by: 10) // 10 sequence chunks
Channel.fromPath('misc/sample.fa')
       .splitFasta(record: [id: true, seqString: true])
       .filter { record -> record.id =~ /^ENST0.*/ } // Record objects
       .subscribe { record -> println record.seqString }
```

* By default chunks are kept in memory.
* Cache to file using `file: true`.
* Many parameters available.

`splitFastq`:

```groovy
Channel.fromPath('misc/sample.fastq')
       .splitFastq(by: 10) // 10 sequence chunks
Channel.fromPath('misc/sample.fastq')
       .splitFastq(record: true) // Record objects
       .println { record -> record.readHeader }
```

* Split paired-end read pairs:

```groovy
Channel.fromFilePairs('/my/data/SRR*_{1,2}.fastq', flat:true)
       .splitFastq(by: 100_000, pe:true, file:true)
       .println()
```

* By default chunks are kept in memory.
* Cache to file using `file: true`.
* Many parameters available.

`splitText`:

* Split multi-line strings or text file items.

```groovy
Channel.fromPath('/some/path/*.txt').splitText() // Single lines.
Channel.fromPath('/some/path/*.txt').splitText(by: 10) // 10 line chunks.
Channel.fromPath('/some/path/*.txt').splitText(by: 10) { it.toUpperCase() }
```

### Combining operators

`join`:

* Join together items where first occurrence of each key matches.

```groovy
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right).println()
```
```console
[Z, 3, 6]
[Y, 2, 5]
[X, 1, 4]
```
```groovy
left = Channel.from(['X', 1], ['Y', 2], ['Z', 3], ['P', 7])
right = Channel.from(['Z', 6], ['Y', 5], ['X', 4])
left.join(right, remainder: true).println()
```
```console
[Y, 2, 5]
[Z, 3, 6]
[X, 1, 4]
[P, 7, null]
```

`merge`:

```groovy
odds  = Channel.from([1, 3, 5, 7, 9]);
evens = Channel.from([2, 4, 6]);
odds.merge(evens).println()
```
```console
[1, 2]
[3, 4]
[5, 6]
```
```groovy
odds.merge(evens) { a, b -> tuple(b*b, a) }.println()
```

`mix`:

```groovy
c1 = Channel.from(1, 2, 3)
c2 = Channel.from('a', 'b')
c3 = Channel.from('z')
c1.mix(c2,c3)
  .subscribe onNext: { println it }, onComplete: { println 'Done' } // 1 2 3 'a' 'b' 'z' Done
```

* Order cannot be guaranteed.

`cross`:

```groovy
source = Channel.from([1, 'alpha'], [2, 'beta'])
target = Channel.from([1, 'x'], [1, 'y'], [1, 'z'], [2, 'p'], [2, 'q'], [2, 't'])
source.cross(target).subscribe { println it }
```
```console
[ [1, alpha], [1, x] ]
[ [1, alpha], [1, y] ]
[ [1, alpha], [1, z] ]
[ [2, beta],  [2, p] ]
[ [2, beta],  [2, q] ]
[ [2, beta],  [2, t] ]
```

`collectFile`:

* Save items to file and output file.

```groovy
Channel.from('alpha', 'beta', 'gamma')
       .collectFile(name: 'sample.txt', newLine: true)
       .subscribe {
           println "Entries are saved to file: $it"
           println "File content is: ${it.text}"
       }
```
```groovy
Channel.from('Hola', 'Ciao', 'Hello', 'Bonjour', 'Halo')
      .collectFile() { item ->
          [ "${item[0]}.txt", item + '\n' ]
      }
      .subscribe {
          println "File ${it.name} contains:"
          println it.text
      }
```
```console
File 'B.txt' contains:
Bonjour

File 'C.txt' contains:
Ciao

File 'H.txt' contains:
Halo
Hola
Hello
```
```groovy
Channel.from('Z'..'A')
       .collectFile(name:'result', sort: true, newLine: true)
       .subscribe { println it.text } // A B C ... Z
```

* Collect and sort all sequences in a FASTA file, shortest to longest:

```groovy
Channel
     .fromPath('/data/sequences.fa')
     .splitFasta(record: [id: true, sequence: true])
     .collectFile(name:'result.fa', sort: { it.size() }) { it.sequence }
     .subscribe { println it.text }
```

`combine` (cartesian product):

```groovy
numbers = Channel.from(1, 2, 3)
words = Channel.from('hello', 'ciao')
numbers.combine(words).println()
```
```console
[1, hello]
[2, hello]
[3, hello]
[1, ciao]
[2, ciao]
[3, ciao]
```
```groovy
left = Channel.from(['A', 1], ['B', 2], ['A', 3])
right = Channel.from(['B','x'], ['B','y'], ['A','z'], ['A', 'w'])
left.combine(right, by: 0).println()
```
```console
[A, 1, z]
[A, 3, z]
[A, 1, w]
[A, 3, w]
[B, 2, x]
[B, 2, y]
```

`concat`:

* Concatenate in known order, unlike `mix`.

```groovy
a = Channel.from('a', 'b', 'c')
b = Channel.from(1, 2, 3)
c = Channel.from('p', 'q')
c.concat(b, a).subscribe { println it } // p q 1 2 3 a b c
```

### Forking operators

`branch`:

* Experimental.

```groovy
Channel
    .from(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it > 10
    }
    .set { result }
result.small.view { "$it is small" }
result.large.view { "$it is large" }
```
```console
1 is small
2 is small
3 is small
40 is large
50 is large
```
```groovy
Channel
    .from(1, 2, 3, 40, 50)
    .branch {
        small: it < 10
        large: it < 50
        other: true
    }
```
```groovy
Channel
    .from(1, 2, 3, 40, 50)
    .branch {
        foo: it < 10
            return it + 2
        bar: it < 50
            return it - 2
        other: true
            return 0
}
```

* `branchCriteria` built-in method.

```groovy
def criteria = branchCriteria {
                small: it < 10
                large: it > 10
               }
Channel.from(1, 2, 30).branch(criteria).set { ch1 }
Channel.from(10, 20, 1).branch(criteria).set { ch2 }
```

`multiMap`:

* Experimental.

```groovy
Channel.from(1, 2, 3)
       .multiMap { it ->
           foo: it + 1
           bar: it * it
           }
       .set { result }
result.foo.view { "foo $it" }
result.bar.view { "bar $it" }
```
```console
foo 2
foo 3
foo 4
bar 1
bar 4
bar 9
```

* Shorthand if both channels receive same source items.

```groovy
Channel.from(1,2,3)
       .multiMap { it -> foo: bar: it }
       .set { result }
```
```groovy
def criteria = multiMapCriteria {
                  small: it < 10
                  large: it > 10
                }
Channel.from(1,2,30).multiMap(criteria).set { ch1 }
Channel.from(10,20,1).multiMap(criteria).set { ch2 }
```

* `multiMapCriteria` built-in method.

`into`:

```groovy
Channel.from('a', 'b', 'c')
       .into{ foo; bar }

foo.println{ "Foo emit: " + it }
bar.println{ "Bar emit: " + it }
```
```console
Foo emit: a
Foo emit: b
Foo emit: c
Bar emit: a
Bar emit: b
Bar emit: c
```
```groovy
(foo, bar) = Channel.from( 'a','b','c').into(2)
foo.println{ "Foo emit: " + it }
bar.println{ "Bar emit: " + it }
```

`tap`:

```groovy
log1 = Channel.create().subscribe { println "Log 1: $it" }
log2 = Channel.create().subscribe { println "Log 2: $it" }
Channel.from('a', 'b', 'c')
       .tap(log1)
       .map { it * 2 }
       .tap(log2)
       .subscribe { println "Result: $it" }
```
```console
Log 1: a
Log 1: b
Log 1: c

Log 2: aa
Log 2: bb
Log 2: cc

Result: aa
Result: bb
Result: cc
```
```groovy
Channel.from('a', 'b', 'c')
       .tap { log1 }
       .map { it * 2 }
       .tap { log2 }
       .subscribe { println "Result: $it" }
log1.subscribe { println "Log 1: $it" }
log2.subscribe { println "Log 2: $it" }
```

### Maths operators

`count`:

```groovy
Channel.from(9, 1, 7, 5).count() // 4
Channel.from(4, 1, 7, 1, 1).count(1) // 3
Channel.from('a','c','c','q','b').count (~/c/) // 2
Channel.from('a', 'c', 'c', 'q', 'b').count { it <= 'c' } // 4
```

`countBy`:

```groovy
Channel.from('x', 'y', 'x', 'x', 'z', 'y').countBy() // [x:3, y:2, z:1]
Channel.from('hola', 'hello', 'ciao', 'bonjour', 'halo').countBy { it[0] } // [h:3, c:1, b:1]
```

`min`:

```groovy
Channel.from(8, 6, 2, 5).min() // 2
Channel.from("hello", "hi", "hey").min { it.size() } // "hi"
Channel.from("hello", "hi", "hey").min { a,b -> a.size() <=> b.size() } // "hi"
```

* `A <=> B` return -1 if A < B, 0 if A == B, 1 if A > B.

`max`:

```groovy
Channel.from(8, 6, 2, 5).max() // 8
Channel.from("hello", "hi", "hey").max { it.size() } // "hello"
Channel.from("hello", "hi", "hey").max { a,b -> a.size() <=> b.size() } // "hello"
```

`sum`:

```groovy
Channel.from(8, 6, 2, 5).sum() // 21
Channel.from(4, 1, 7, 5).sum { it * it } // 91
```

`toInteger`:

```groovy
Channel.from( '1', '7', '12' ).toInteger().sum()
```

### Other operators

`dump`:

* Print items emitted if `-dump-channels <CHANNEL>` is specified on command-line.
* Debugging.

```groovy
Channel.from(1, 2, 3).map { it+1 }.dump(tag: 'foo')
Channel.from(1, 2, 3).map { it^2 }.dump(tag: 'bar')
```

`set`:

* Assign channel to variable.
* Can only be the last operator in a chain.

```groovy
Channel.from(10, 20, 30).set { my_channel }
```

* Equivalent to:

```groovy
my_channel = Channel.from(10, 20, 30)
```

`ifEmpty`:

* Default value if channel is empty:

```groov
Channel.from(1, 2, 3).ifEmpty('Hello').println() // 1 2 3
Channel.empty().ifEmpty('Hello') .println() // Hello
```

`print`:

* Print to standard output:
* Can only be the last operator in a chain.

```groovy
Channel.from('foo', 'bar', 'baz', 'qux').print { it.toUpperCase() + ' ' }
```

`println`:

* Can only be the last operator in a chain.
* As for `print` but with a newline added to each item.

```groovy
Channel.from('foo', 'bar', 'baz', 'qux').println { "~ $it" }
```

`view`:

* Print to *console* standard output with a newline added to each item.
* Returns a newly-created channel.
* Unlinke `print`, `view` can be chained like other operators.

```groovy
Channel.from(1, 2, 3).view() // 1 2 3
Channel.from(1, 2, 3).map { it -> [it, it*it] }.view { num, sqr -> "Square of: $num is $sqr" }
```

`close`:

* Close channel and cause downstream processes or operators to stop.
* Not usually used.

---

## Executors

https://www.nextflow.io/docs/latest/executor.html

### Local

* Run in local host.
* Parallelise using multiple threads and multi-core architecture on CPU.

---

## Configuration

https://www.nextflow.io/docs/latest/config.html

### Configuration file

`nextflow.config`:

* Current directory.
* Script base directory, if not the current directory.
* `$HOME/.nextflow/config`

Custom configuration files:

* `-c <config file>` takes precedence over the above.
* `-C <config file>` - only use the given file and ignore the others.

Plain-text file:

* Numbers, quoted strings, `true`, `false`

```
name = value
```

* Dereference configuration variables.
* Dereference environment variables.

```
propertyOne = 'world'
anotherProp = "Hello $propertyOne"
customPath = "$PATH:/my/app/folder"
```

```
// Comments
/*
 ...
 */
```

```
includeConfig 'path/foo.config'
```

### Scopes

```
alpha.x  = 1
alpha.y  = 'string value..'
beta {
     p = 2
     q = 'another string ..'
}
```

`env`:

* Exported to workflow task execution environment.

```
env.ALPHA = 'some value'
env.BETA = "$HOME/some/path"

env {
     DELTA = 'one more'
     GAMMA = "/my/path:$PATH"
}
```

`params`:

* Accessible in pipeline script.

```
params.custom_param = 123
params.another_param = 'string value .. '

params {
   alpha_1 = true
   beta_2 = 'another string ..'
}
```

`process`:

* Default configuration for processes.
* Use in `process` directive and `executor` section.
* Define multiple platform-specific process configurations.

Many other scopes available e.g. `process`, `executor`, `docker`, `conda`, etc.

### Profiles

Configuration profiles

```
profiles {
    standard {
        process.executor = 'local'
    }
    cluster {
        process.executor = 'sge'
        process.queue = 'long'
        process.memory = '10GB'
    }
    cloud {
        process.executor = 'cirrus'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }
}
```
```console
$ nextflow run <your script> -profile standard, cloud
```

### Environment variables

Environment variables controlling configuration of the Nextflow runtime and the JVM.

---

## Tracing & visualisation

HTML report with summary, resource usage, metrics tasks:

```console
$ nextflow run <pipeline name> -with-report [file name]
```

Process timeline:

```console
$ nextflow run <pipeline name> -with-timeline [file name]
```

Pipeline DAG:

```
$ nextflow run <script-name> -with-dag flowchart.dot|html|pdf|png|svg|gexf
$ nextflow run <script-name> -with-dag flowchart.png
```

---

## Pipeline sharing

BitBucket, GitHub, GitLab.

```console
$ nextflow run foo/bar # GitHub
$ nextflow run foo/bar -hub bitbucket
$ nextflow run foo/bar -hub gitlab
$ nextflow run http://github.com/foo/bar
...
```
```console
$ nextflow run nextflow-io/hello
```

This has a `main.nf` file.

```console
$ nextflow run nextflow-io/hello -r mybranch
$ nextflow run nextflow-io/hello -r v1.1
```

A number of other functions related to revision control and repository hosting are provided.
