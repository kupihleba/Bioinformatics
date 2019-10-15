package main

import (
	. "Bioinformatics/Sequence_alignment/algoritm"
	. "Bioinformatics/Sequence_alignment/utils"
	"bufio"
	"flag"
	"fmt"
	"github.com/pkg/errors"
	"io"
	"os"
	"strconv"
	"strings"
)

func readFile(path string) (string, string) {
	file, err := os.OpenFile(path, os.O_RDONLY, os.ModePerm)
	defer file.Close()
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
	defer file.Close()
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

func main() {
	// Command line arguments
	gapPtr := flag.Int("g", -2, "gap penalty as int")
	inpPtr := flag.String("i", "",
		"input file, containing 2 sequences, separated with a newline")
	outpPtr := flag.String("o", "",
		"output file, write 2 aligned sequences, separated with a newline")
	typePtr := flag.String("t", "default",
		"type of the weight matrix. Possible types DNAFull, BLOSUM62, DEFAULT")
	flag.Parse()

	var (
		seq1, seq2 string
		err        error
	)
	inpFile := strings.TrimSpace(*inpPtr)
	if inpFile != "" {
		seq1, seq2 = readFile(inpFile)
		check(err)
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

	alignedSeq1, alignedSeq2, score, err := engine.NeedlemanWunsch(seq1, seq2)
	check(err)

	outpFile := strings.TrimSpace(*outpPtr)
	if outpFile != "" {
		writeSeqToFile(outpFile, alignedSeq1, alignedSeq2, score)
	} else {
		fmt.Printf("Aligned seq1:\t%s\nAligned seq2:\t%s\nScore: %d",
			Prettify(alignedSeq1, 100), Prettify(alignedSeq2, 100), score)
	}
}
