package main

import (
	. "Bioinformatics/Sequence_alignment/algorithm"
	. "Bioinformatics/Sequence_alignment/utils"
	"bufio"
	"flag"
	"fmt"
	"github.com/pkg/errors"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

func readTemplate(path string) string {
	file, err := os.OpenFile(path, os.O_RDONLY, os.ModePerm)
	defer func() {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}()

	check(err)
	reader := bufio.NewReader(file)

	seq1, err := reader.ReadString('\n')
	if err != io.EOF {
		check(err)
	}

	return strings.TrimSpace(seq1)
}

func readFile(path string) (string, string) {
	file, err := os.OpenFile(path, os.O_RDONLY, os.ModePerm)
	defer func() {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}()

	check(err)
	reader := bufio.NewReader(file)

	seq1, err := reader.ReadString('\n')
	check(err)
	seq2, err := reader.ReadString('\n')
	if err != io.EOF {
		check(err)
	}

	return strings.TrimSpace(seq1), strings.TrimSpace(seq2)
}

func writeSeqToFile(path string, seq1 string, seq2 string, score int) {
	file, err := os.Create(path)
	defer func() {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}()
	check(err)
	_, err = file.WriteString(seq1)
	check(err)
	_, err = file.WriteString("\n")
	check(err)
	_, err = file.WriteString(seq2)
	check(err)
	_, err = file.WriteString("\n")
	check(err)
	_, err = file.WriteString(strconv.Itoa(score))
	check(err)
}

func getScoreFuncAndPenalty(funcType string) (ScoreFuncType, int, error) {
	funcType = strings.ToLower(strings.TrimSpace(funcType))
	switch funcType {
	case "default":
		return ScoreDefault, -2, nil
	case "blosum62":
		return ScoreBLOSUM62, -4, nil
	case "dnafull":
		return ScoreDNAFull, -4, nil
	default:
		return nil, -2, errors.New("Unknown matrix type")
	}
}

func check(err error) {
	if err != nil {
		panic(err)
	}
}

func isFlagPassed(name string) bool {
	found := false
	flag.Visit(func(f *flag.Flag) {
		if f.Name == name {
			found = true
		}
	})
	return found
}

func readFastaFilePart(reader *bufio.Reader, maxSize int) ([]string, bool) {
	arr := make([]string, maxSize)
	i := 0
	var sb strings.Builder
	var str string
	var err error
	for err == nil && i < maxSize {
		str, err = reader.ReadString('\n')
		str = strings.TrimSpace(str)
		prefix := strings.HasPrefix(str, ">")
		if prefix && sb.Len() > 0 {
			arr[i] = sb.String()
			i++
			sb.Reset()
		} else if !prefix {
			sb.WriteString(str)
		}
	}
	if sb.Len() > 0 {
		arr[i] = sb.String()
		i++
	}
	return arr[:i], err == io.EOF
}

func goFastaCompute(template string, sequences []string, engine AlignEngine) DataChunk {
	const PART_SIZE = 1_000
	if len(sequences) < PART_SIZE {
		if len(sequences) == 0 {
			return DataChunk{template, "", math.MinInt64}
		}
		_, _, score, index := engine.MultiAlignSequences(template, sequences)
		return DataChunk{
			str1:  template,
			str2:  sequences[index],
			score: score,
		}
	}
	workers := 0
	ch := make(chan DataChunk)
	i := 0
	j := i + PART_SIZE
	for i < len(sequences) {
		if j >= len(sequences) {
			j = len(sequences) - 1
		}
		workers++
		scope := sequences[i:j]

		go func(template string, scope []string) {
			if len(scope) > 0 {
				_, _, score, index := engine.MultiAlignSequences(template, scope)
				ch <- DataChunk{
					str1:  template,
					str2:  scope[index],
					score: score,
				}
			} else {
				ch <- DataChunk{"", "", math.MinInt64}
			}
		}(template, scope)
		if i == j {
			break
		}
		i = j
		j = i + PART_SIZE
	}

	var resSequences []string
	for i := 0; i < workers; i++ {
		data := <-ch
		if len(data.str2) > 0 {
			resSequences = append(resSequences, data.str2)
		}
		//fmt.Printf("str1:\t%s\nstr2:\t%s\nscore:\t%d\n", data.str1, template, data.score)
	}
	return goFastaCompute(template, resSequences, engine)
}

func goFasta(path string, template string, engine AlignEngine) (string, string, int, int) {
	file, err := os.OpenFile(path, os.O_RDONLY, os.ModePerm)
	defer func() {
		err := file.Close()
		if err != nil {
			log.Fatal(err)
		}
	}()
	check(err)
	reader := bufio.NewReader(file)
	isEOF := false
	var sequences []string
	const PART_SIZE = 100_000
	workers := 0
	var resSequences []string
	for !isEOF { // So, let's read file by parts and spawn workers for each part
		sequences, isEOF = readFastaFilePart(reader, PART_SIZE)
		fmt.Printf("Computation stage %d\n", workers)
		workers++
		seqCandidate := goFastaCompute(template, sequences, engine).str2
		if len(seqCandidate) > 0 {
			resSequences = append(resSequences, seqCandidate)
		}
	}

	//fmt.Printf("Spawned %d workers\n", workers)
	//for i := 0; i < workers; i++ {
	//	data := <-ch
	//	resSequences = append(resSequences, data.str2)
	//	fmt.Printf("str1:\t%s\nstr2:\t%s\nscore:\t%d\n", data.str1, data.str2, data.score)
	//}
	fmt.Printf("res seq: %v\n", resSequences)
	return engine.MultiAlignSequences(template, resSequences)
}

type DataChunk struct {
	str1, str2 string
	score      int
}

func main() {
	// Command line arguments
	gapPtr := flag.Int("g", -2, "gap penalty as int")
	inpPtr := flag.String("i", "",
		"input file, containing 2 sequences, separated with a newline")
	outpPtr := flag.String("o", "",
		"output file, write 2 aligned sequences, separated with a newline")
	typePtr := flag.String("t", "default",
		"type of the weight matrix. Possible types DNAFull, BLOSUM62, DEFAULT")
	algoPtr := flag.String("algo", "Needleman-Wunsch",
		"Chose the alignment algorithm (Needleman-Wunsch|Smith-Waterman|Hirschberg|FASTA)")
	//multiAlignPtr := flag.String("fasta", "",
	//	"Read file in FASTA format and go FASTA!")
	templatePtr := flag.String("templ", "",
		"Template for FASTA alignment")
	flag.Parse()

	algo := strings.TrimSpace(*algoPtr)
	algo = strings.Replace(algo, "-", "", -1)
	algo = strings.ToLower(algo)

	var (
		seq1, seq2 string
		err        error
	)
	inpFile := strings.TrimSpace(*inpPtr)
	if inpFile != "" && algo != "fasta" {
		seq1, seq2 = readFile(inpFile)
		check(err)
	} else if *templatePtr != "" {
	} else {
		_, _ = fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
		flag.PrintDefaults()
		return
	}
	scoreFunc, gapPenalty, err := getScoreFuncAndPenalty(*typePtr)
	check(err)

	if isFlagPassed("g") {
		gapPenalty = *gapPtr
	}

	// Engine setup
	engine := NewAlignEngine(scoreFunc, gapPenalty)

	seq1 = strings.ToUpper(seq1)
	seq2 = strings.ToUpper(seq2)

	var (
		alignedSeq1, alignedSeq2 string
		score                    int
	)

	switch algo {
	case "hirschberg":
		alignedSeq1, alignedSeq2 = engine.Hirschberg(seq1, seq2)
		break
	case "smithwaterman":
		alignedSeq1, alignedSeq2, score = engine.SmithWaterman(seq1, seq2)
		break
	case "needlemanwunsch":
		alignedSeq1, alignedSeq2, score = engine.NeedlemanWunsch(seq1, seq2)
		break
	case "fasta":
		if *templatePtr == "" {
			panic("Pass FASTA template!")
		}
		template := readTemplate(*templatePtr)
		if inpFile == "" {
			panic("No input file specified!")
		}
		alignedSeq1, alignedSeq2, score, _ = goFasta(inpFile, template, engine)
		break
	default:
		panic("Unknown algorithm! Available options = Needleman-Wunsch | Smith-Waterman | Hirschberg | FASTA")
	}
	check(err)

	outpFile := strings.TrimSpace(*outpPtr)
	if outpFile != "" {
		writeSeqToFile(outpFile, alignedSeq1, alignedSeq2, score)
	} else {
		fmt.Printf("Aligned seq1:\t%s\nAligned seq2:\t%s\nScore: %d",
			Prettify(alignedSeq1, 100), Prettify(alignedSeq2, 100), score)
	}
}

//Aligned seq1:   LKMYGIVTTVKLANKMIQNEKFEVWDILDEVIHEHPIL
//Aligned seq2:   LKMYGIVPTVKLANKMIQNEKPEVWDILDEVIHEHPIL
//Score: 34
