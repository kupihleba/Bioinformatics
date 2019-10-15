// Helpful methods that are not in std library
package utils

import "strings"

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

func Max(numbers ... int) int {
	max_numb := numbers[0]
	for i := 1; i < len(numbers); i++ {
		if max_numb < numbers[i] {
			max_numb = numbers[i]
		}
	}
	return max_numb
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
