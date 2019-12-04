#ifndef WIGGLETOSITE_HPP
#define WIGGLETOSITE_HPP

#include "../../../GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/includes/gmml.hpp"

class WiggleToSite
{
public:

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////

    WiggleToSite();
    WiggleToSite(MolecularModeling::Assembly &moving_assembly, MolecularModeling::Assembly &receptor_assembly, MolecularModeling::Assembly &target_assembly, const std::string inputFile);

    //////////////////////////////////////////////////////////
    //                       TYPEDEFS                       //
    //////////////////////////////////////////////////////////

    typedef std::vector< std::pair <Residue *, Residue *> > ResiduePairVector;

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    AtomVector GetMovingAssembly();
    AtomVector GetReceptorAssembly();
    AtomVector GetSortedMovingAtoms();
    AtomVector GetSortedTargetAtoms();

    //////////////////////////////////////////////////////////
    //                       MUTATOR                        //
    //////////////////////////////////////////////////////////

    void SetMovingAssembly(MolecularModeling::Assembly &moving_assembly);
    void SetReceptorAssembly(MolecularModeling::Assembly &receptor_assembly);
    void SetTargetAssembly(MolecularModeling::Assembly &target_assembly);
    //void AddmovingTargetResiduePairs(Residue* moving, Residue* target);

    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////


    //    void WiggleResolveOverlap(double CheckValueFunction(), int interval = 5, int iterations = 5, double acceptableValue = 0.1);
//    void WiggleDistanceOverlap(int interval = 5, int iterations = 5, double acceptableDistance = 1.0);
//    void WigglePermutatorDistanceOverlap(int interval = 5, int iterations = 1, double acceptableDistance = 1.0);
    double CalculateTargetDistance();
    double CalculateOverlaps();
    void WritePdbFile(std::string fileName = "wigglerFile.pdb");
//    void generateRotatableDihedralPermutations(double &lowest_distance, RotatableDihedralVector::iterator rotatableDihedral, RotatableDihedralVector::iterator end, int interval = 5);
    void StashCoordinates();
    void SetStashedCoordinates();
    void WiggleSimpleDistance(int interval = 5, int iterations = 1);
    void WigglePermutatorDistance(int interval = 5, int iterations = 50, int dihedralsPerIteration = 3);
    //////////////////////////////////////////////////////////
    //                       DISPLAY FUNCTION               //
    //////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //                       OPERATORS                      //
    //////////////////////////////////////////////////////////

private:

    //////////////////////////////////////////////////////////
    //                    PRIVATE FUNCTIONS                 //
    //////////////////////////////////////////////////////////
    void GenerateRotatableDihedralPermutationsDistance(double &lowest_distance, double &lowest_overlap,RotatableDihedralVector::iterator rotatableDihedral, RotatableDihedralVector::iterator end, int interval);

    void ReadInputFile(const std::string inputFile);
    //void GenerateMovingTargetResiduePairs();
    void SortMovingTargetByAtomNameIntoAtomVectors();
    void SortMovingTargetSuperimpositionByAtomNameIntoAtomVectors();
    void SuperimposeToTargetResidues(AtomVector supers, AtomVector targets);
    int generatuniqueID();

    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    MolecularModeling::Assembly* moving_assembly_;
    MolecularModeling::Assembly* receptor_assembly_;
    MolecularModeling::Assembly* target_assembly_;
    ResidueVector movingResidues_;
    ResidueVector targetResidues_;
    ResidueVector superimpositionResidues_;
    ResidueVector superimpositionTargetResidues_;
    ResidueLinkageVector adjustableLinkages_;
    AtomVector sortedMovingAtoms_;
    AtomVector sortedTargetAtoms_;
    AtomVector sortedSuperimpositionAtoms_;
    AtomVector sortedSuperimpositionTargetAtoms_;

};

//std::ostream& operator<<(std::ostream& os, const WiggleToSite&);

#endif // WIGGLETOSITE_HPP
