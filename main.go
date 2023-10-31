package main

import (
	"fmt"
	"os"
	"strings"

	"github.com/Open-Science-Global/poly/io/fasta"
	"github.com/Open-Science-Global/poly/transform/codon"
)

func GenerateCodonTable(file string) codon.Table {
	// input: a list of CDSs from the target organism
	// ouput:

	fmt.Printf("Reading file %s...\n", file)
	cdsSequences := fasta.Read(file)

	// Create a single big string with all the CDSs
	var allCdssFromFile strings.Builder
	for _, cds := range cdsSequences {
		allCdssFromFile.WriteString(cds.Sequence)
	}
	codingRegions := allCdssFromFile.String()

	fmt.Printf("Creating table for %s...\n", file)
	codonTable := codon.GetCodonTable(11)

	fmt.Printf("Optimizing table for %s...\n", file)
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	fmt.Println("Table created and optimized!")
	fmt.Printf("\n")
	return optimizationTable

}

func CodonOptimization(enzymeSequence string, codonTable codon.Table) string {
	// input

	// Optimize sequence using the protein sequence and codon table
	optimizedSequence, _ := codon.Optimize(enzymeSequence, codonTable)

	// Lets check if the codon optimization actually works by making some checks:
	// First one is if both codon sequences are different
	if optimizedSequence == enzymeSequence {
		fmt.Println("Both sequences are equal, some problem occur. They should be different because one is optimized. Checks what happened and run again.")
		os.Exit(0)
	}

	// Check if both translated sequences are equal
	protein, _ := codon.Translate(optimizedSequence, codon.GetCodonTable(11))
	if protein != enzymeSequence {
		fmt.Println("These protein sequences aren't equal, some problem occur. They should be equal because codon optimization don't change any aminoacid.")
		os.Exit(0)
	}
	return optimizedSequence
}

func exportSequencesAsFasta(sequences []string, enzymes []fasta.Fasta, outputFilename string) {
	var fastas []fasta.Fasta
	for index, sequence := range sequences {
		data := fasta.Fasta{enzymes[index].Name, sequence}
		fmt.Println(data)
		fastas = append(fastas, data)
	}
	fasta.Write(fastas, outputFilename)
}

func main() {
	codonTable := GenerateCodonTable("PichiaPastoris-Genome.fasta")
	enzymes := fasta.Read("enzymes.fasta")
	fmt.Print(len(enzymes))
	fmt.Println(enzymes[0])
	fmt.Println(enzymes[1])

	var enzymesCodonOptimized []string
	for _, enzyme := range enzymes {
		enzymesCodonOptimized = append(enzymesCodonOptimized, CodonOptimization(enzyme.Sequence, codonTable))
	}

	fmt.Print(len(enzymesCodonOptimized))

	for _, sequence := range enzymesCodonOptimized {
		fmt.Println(sequence)
	}
	// run this to set the function to export your optimized sequence in a fasta file
	// input: the output from CodonOptimization()
	// output: one fasta file with all optimized sequences

	exportSequencesAsFasta(enzymesCodonOptimized, enzymes, "OptmizedSequencesLinker.fasta")
}
