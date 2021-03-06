// Helpful methods that are not in std library
package utils

import (
	"strings"
)

func ReverseStr(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func ReverseArray(arr []int) {
	l := len(arr)
	for i := 0; i < l/2; i++ {
		arr[i], arr[l-i-1] = arr[l-i-1], arr[i]
	}
}

// Function returns 2 values: first - max value, second - index of max value
func Max(numbers ...int) (int, int) {
	max_numb := numbers[0]
	pos := 0
	for i := 1; i < len(numbers); i++ {
		if max_numb < numbers[i] {
			max_numb = numbers[i]
			pos = i
		}
	}
	return max_numb, pos
}

func Min(numbers ...int) (int, int) {
	min_numb := numbers[0]
	pos := 0
	for i := 1; i < len(numbers); i++ {
		if min_numb > numbers[i] {
			min_numb = numbers[i]
			pos = i
		}
	}
	return min_numb, pos
}

// Iterates over arrays, in search of the max sum of the elements with the same indices
func SumAndMax(arr1 []int, arr2 []int) int {
	if len(arr1) != len(arr2) {
		panic("SumAndMax: arrays must be of the same length!")
	} else if len(arr1) == 0 {
		panic("SumAndMax: arrays must not be empty!")
	}
	iMax := 0
	valMax := arr1[0] + arr2[0]
	for i := 1; i < len(arr1); i++ {
		maxCandidate := arr1[i] + arr2[i]
		if valMax < maxCandidate {
			iMax = i
			valMax = maxCandidate
		}
	}
	return iMax
}

// Inserts \n to fit lines into provided length
func Prettify(sequence string, length int) string {
	var sb strings.Builder
	for len(sequence) > length {
		sb.WriteString(sequence[:length])
		sb.WriteByte('\n')
		sequence = sequence[length:]
	}
	sb.WriteString(sequence)
	return sb.String()
}
