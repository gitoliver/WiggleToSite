#include <iostream>
#include "../../../gems/gmml/includes/gmml.hpp"
#include "../includes/io.h"
#include "../includes/wiggleTosite.hpp"

//double CalculateAverageDistanceBetweenAtomVectors(AtomVector A, AtomVector B);
int main(int argc, char *argv[])
{
    std::cout << "Hello Cruel World!" << std::endl;

    std::string workingDirectory = Find_Program_Working_Directory(); // Default behaviour.
    std::string installationDirectory = Find_Program_Installation_Directory();

    std::string inputFile, movingPDB, overlapPDB, targetPDB;
    if ( (argc == 4) || (argc == 5) )
    {
        inputFile = workingDirectory + "/" + argv[1];
        movingPDB = argv[2];
        overlapPDB = argv[3];
        targetPDB = argv[3];
        if (argc == 5)
            targetPDB = argv[4];
    }
    else
    {
        std::cerr << "There ain't enough arguments (argc is " << argc << ")\n";
        std::cerr << "Usage: " << argv[0] << " inputFile movingPDB overlapPDB [targetPDB]. Note: A separate targetPDB is optional, target residues can be in overlapPDB.\n";
        std::exit (1);
    }

    std::cout << "Working directory is " << workingDirectory << "\n";
    std::cout << "Install directory is " << installationDirectory << "\n"; // Assumes folder with glycans inside is present.

    MolecularModeling::Assembly movingAss(movingPDB, gmml::PDB);
    MolecularModeling::Assembly overlapAss(overlapPDB, gmml::PDB);
    MolecularModeling::Assembly targetAss(targetPDB, gmml::PDB);
    movingAss.BuildStructureByDistance(4, 1.6); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded. Nope that did bad things
    overlapAss.BuildStructureByDistance(4, 1.6); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded. Nope that did bad things
    targetAss.BuildStructureByDistance(4, 1.6); // 4 threads, 1.91 cutoff to allow C-S in Cys and Met to be bonded. Nope that did bad things


    WiggleToSite wiggler(movingAss, overlapAss, targetAss, inputFile);

    //wiggler.WigglePermutatorDistance();
  //  wiggler.WiggleSimpleDistance();
    wiggler.WigglePermutatorDistance(1, 30, 3);
//    wiggler.WiggleDistanceOverlap(5, 5, 1.0);
//    wiggler.WiggleDistanceOverlap(2, 2, 1.0);
//    wiggler.WiggleDistanceOverlap(1, 1, 0.5);


    PdbFileSpace::PdbFile *outputPdbFile = movingAss.BuildPdbFileStructureFromAssembly(-1,0);
    outputPdbFile->Write(workingDirectory + "/wiggled.pdb");

    // Add beads. They make the overlap calculation faster.
    //    std::cout << "Add_Beads"  << std::endl;
    //    beads::Add_Beads(glycoprotein, glycosites);
    //    glycoprotein_builder::UpdateAtomsThatMoveInLinkages(glycosites); // Must update to include the beads

    return 0;
}

//double CalculateAveragePerAtomDistanceBetweenSortedAtomVectors(AtomVector atomsA, AtomVector atomsB)
//{
//    if (atomsA.size() != atomsB.size() || atomsA.empty())
//    {
//        std::cerr << "Problem in CalculateAveragePerAtomDistanceBetweenSortedAtomVectors. Exiting" << std::endl;
//        std::exit(1);
//    }
//    double distanceSum = 0.0;
//    for (int i = 0; i < atomsA.size(); ++i)
//    {
//        distanceSum += atomsA.at(i)->GetDistanceToAtom(atomsB.at(i));
//    }
//    return (distanceSum / atomsA.size());
//}


