#include "../includes/wiggleTosite.hpp"
#include "../includes/io.h" // split function

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
/////////////////////////////////////////////////////////

WiggleToSite::WiggleToSite()
{

}

WiggleToSite::WiggleToSite(MolecularModeling::Assembly &moving_assembly, MolecularModeling::Assembly &receptor_assembly, MolecularModeling::Assembly &target_assembly, const std::string inputFile)
{
    this->SetMovingAssembly(moving_assembly);
    this->SetReceptorAssembly(receptor_assembly);
    this->SetTargetAssembly(target_assembly);
    this->ReadInputFile(inputFile);
    this->SortMovingTargetByAtomNameIntoAtomVectors();
}

//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////

AtomVector WiggleToSite::GetMovingAssembly()
{
    return moving_assembly_->GetAllAtomsOfAssembly();
}

AtomVector WiggleToSite::GetReceptorAssembly()
{
    return receptor_assembly_->GetAllAtomsOfAssembly();
}

AtomVector WiggleToSite::GetSortedMovingAtoms()
{
    return sortedMovingAtoms_;
}

AtomVector WiggleToSite::GetSortedTargetAtoms()
{
    return sortedTargetAtoms_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////

void WiggleToSite::SetMovingAssembly(MolecularModeling::Assembly  &moving_assembly)
{
    moving_assembly_ = &moving_assembly;
}

void WiggleToSite::SetReceptorAssembly(MolecularModeling::Assembly &receptor_assembly)
{
    receptor_assembly_ = &receptor_assembly;
}

void WiggleToSite::SetTargetAssembly(MolecularModeling::Assembly &target_assembly)
{
    target_assembly_ = &target_assembly;
}

//////////////////////////////////////////////////////////
//                       FUNCTIONS                      //
//////////////////////////////////////////////////////////

double WiggleToSite::CalculateTargetDistance()
{
    AtomVector atomsA = this->GetSortedMovingAtoms();
    AtomVector atomsB = this->GetSortedTargetAtoms();
    double distanceSum = 0.0;
    for (int i = 0; i < atomsA.size(); ++i)
    {
        distanceSum += atomsA.at(i)->GetDistanceToAtom(atomsB.at(i));
    }
    return (distanceSum / atomsA.size());
}

double WiggleToSite::CalculateOverlaps()
{
    // Maybe more efficient to get AtomVectors during construction and have them as attributes.
    return gmml::CalculateAtomicOverlaps(receptor_assembly_->GetAllAtomsOfAssembly(), moving_assembly_->GetAllAtomsOfAssembly());
}

void WiggleToSite::WritePdbFile(std::string fileName)
{
    PdbFileSpace::PdbFile *outputPdbFile = moving_assembly_->BuildPdbFileStructureFromAssembly(-1,0);
    std::stringstream stream;
    stream << std::setfill('0') << std::setw(6) << this->generatuniqueID(); // Gives me 00001, 00002, etc so VMD loads them in order.
    stream << "_" << fileName;
    outputPdbFile->Write(stream.str());
    return;
}

void WiggleToSite::StashCoordinates()
{
    AtomVector atoms = moving_assembly_->GetAllAtomsOfAssembly();
    for(auto &atom : atoms)
    {
        atom->AddCoordinate(new GeometryTopology::Coordinate(atom->GetCoordinate())); // push back the currect coordinates onto the end of the coordinates vector.
    }
    return;
}

void WiggleToSite::SetStashedCoordinates() // Each time a better structure is found, the coords are pushed back. Last one pushed back is the best! (No tacksie backsies).
{
    AtomVector atoms = moving_assembly_->GetAllAtomsOfAssembly();
    for(auto &atom : atoms)
    {
        GeometryTopology::Coordinate *first_coords = atom->GetCoordinate();
        GeometryTopology::Coordinate *last_coords = atom->GetCoordinates().back();
        *first_coords = *last_coords; // Biblical
    }
    return;
}

void WiggleToSite::WiggleSimpleDistance(int interval, int iterations)
{
    for(int i = 1; i <= iterations; ++i)
    {
        if (i > 1)
        { // Give them the old razzle dazzle
            std::random_shuffle(adjustableLinkages_.begin(), adjustableLinkages_.end());
        }
        std::cout << "Iteration: " << i << "\n";
        for(auto &linkage : adjustableLinkages_)
        {
            RotatableDihedralVector rotatable_bond_vector = linkage.GetRotatableDihedrals();
            std::cout << "Simple Working on " << linkage.GetResidues().at(0)->GetId() << "---" << linkage.GetResidues().at(1)->GetId()
                      << " which has " << rotatable_bond_vector.size() << " rotatble bonds." << std::endl;
            double lowest_overlap = this->CalculateOverlaps();
            double lowest_distance = this->CalculateTargetDistance();
            this->GenerateRotatableDihedralPermutationsDistance(lowest_distance, lowest_overlap, rotatable_bond_vector.begin(), rotatable_bond_vector.end(), interval);
        }
    }
    return;
}

void WiggleToSite::WigglePermutatorDistance(int interval, int iterations, int dihedralsPerIteration)
{
    RotatableDihedralVector all_rotatable_bonds;
    all_rotatable_bonds.reserve(adjustableLinkages_.size()*2); // EFFICIENCY, we must have it?
    for(auto &linkage : adjustableLinkages_)
    {
        RotatableDihedralVector current_bonds = linkage.GetRotatableDihedrals();
        all_rotatable_bonds.insert(all_rotatable_bonds.end(), current_bonds.begin(), current_bonds.end());
    }
    for(int i = 1; i <= iterations; ++i)
    {
        if (i > 1)
        { // Give them the old razzle dazzle
            std::random_shuffle(all_rotatable_bonds.begin(), all_rotatable_bonds.end());
        }
        std::cout << "Iteration: " << i << "\n";
        double lowest_overlap = this->CalculateOverlaps();
        double lowest_distance = this->CalculateTargetDistance();
        RotatableDihedralVector::iterator nthElement = (all_rotatable_bonds.begin() + dihedralsPerIteration);
        this->GenerateRotatableDihedralPermutationsDistance(lowest_distance, lowest_overlap, all_rotatable_bonds.begin(), nthElement, interval);
    }
    return;
}

//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

void WiggleToSite::GenerateRotatableDihedralPermutationsDistance(double &lowest_distance, double &lowest_overlap, RotatableDihedralVector::iterator rotatableDihedral, RotatableDihedralVector::iterator end, int interval)
{
    std::vector<double> allAngles = rotatableDihedral->GetAllPossibleAngleValues();
    for(auto &angle : allAngles)
    {
    //    std::cout << angle << ", ";
    //    this->WritePdbFile("all.pdb");
        rotatableDihedral->SetDihedralAngle(angle);
        if(std::next(rotatableDihedral) != end)
        {
            WiggleToSite::GenerateRotatableDihedralPermutationsDistance(lowest_distance, lowest_overlap, std::next(rotatableDihedral), end, interval);
        }
        else // At the end of the vector of RotatableDihedrals, so finished setting a new permutant.
        {
            double current_overlap = this->CalculateOverlaps();
            double current_distance = this->CalculateTargetDistance();
            if ( (current_overlap < lowest_overlap) || ((current_overlap <= lowest_overlap) && ( current_distance < lowest_distance)) )
            {  // Always stash if lower overlap, only stash if distance is lower and overlap is the same or lower
                std::cout << "Stashing coords with\noverlaps: " << current_overlap << "\ndistance: " << current_distance << "\n";
                this->StashCoordinates();
                lowest_overlap = current_overlap;
                lowest_distance = current_distance;
                this->WritePdbFile("best.pdb");
            }
        }
    }
    //std::cout << std::endl;
    this->SetStashedCoordinates(); // The last one stashed will be the best.
    return;
}

//void glycoprotein_builder::Read_Input_File(GlycosylationSiteVector &glycosites, std::string &proteinPDB, std::string &glycanDirectory, const std::string working_Directory)
void WiggleToSite::ReadInputFile(const std::string inputFile)
{
    // Read input file
    std::string buffer;
    std::ifstream infile (inputFile);
    if (!infile)
    {
        std::cerr << "Uh oh, input file could not be opened for reading!" << std::endl;
        std::exit(1);
    }
    while (infile) // While there's still stuff left to read
    {
        std::string strInput;
        getline(infile, strInput);
        // Format should be a new entry on each line, without the quotes: "?_222, ?_223". Where ? indicates no chain ID info. Ends with a line containing "END".
        if(strInput == "AdjustableLinkages:")
        {
            getline(infile, buffer);
            while(buffer != "END")
            {
                StringVector splitLine = split(buffer, ',');
                if(splitLine.size() == 2)
                {
                    adjustableLinkages_.emplace_back(selection::FindResidue(*moving_assembly_, splitLine.at(0)), selection::FindResidue(*moving_assembly_, splitLine.at(1))); // Creates ResidueLinkage instance on the vector. Love it.
                }
                else if(splitLine.size() == 3)
                {
                    if(splitLine.at(2) == "noReverse")
                    {
                        adjustableLinkages_.emplace_back(selection::FindResidue(*moving_assembly_, splitLine.at(0)), selection::FindResidue(*moving_assembly_, splitLine.at(1)), false); // Creates ResidueLinkage instance on the vector. Love it.

                    }
                }
                getline(infile, buffer);
            }
        }
        if(strInput == "MovingResidues:")
        {
            getline(infile, buffer);
            while(buffer != "END")
            {
                movingResidues_.push_back(selection::FindResidue(*moving_assembly_, buffer));
                getline(infile, buffer);
            }
        }
        if(strInput == "TargetResidues:")
        {
            getline(infile, buffer);
            while(buffer != "END")
            {
                targetResidues_.push_back(selection::FindResidue(*target_assembly_, buffer));
                getline(infile, buffer);
            }
        }
    }
    return;
}

void WiggleToSite::SortMovingTargetByAtomNameIntoAtomVectors()
{
    if ( (movingResidues_.size() != targetResidues_.size()) || movingResidues_.empty() )
    {
        std::cerr << "Problem in WiggleToSite::SortMovingTargetByAtomNameIntoAtomVectors()" << std::endl;
        std::exit(1);
    }
    for (int i = 0; i < movingResidues_.size(); ++i)
    {
        AtomVector atomsA = movingResidues_.at(i)->GetAtoms();
        AtomVector atomsB = targetResidues_.at(i)->GetAtoms();
        for (auto &atomA : atomsA)
        {
            for (auto &atomB : atomsB)
            {
                if (atomA->GetName() == atomB->GetName())
                {
                    sortedMovingAtoms_.push_back(atomA);
                    sortedTargetAtoms_.push_back(atomB);
                }
            }
        }
    }
    return;
}

int WiggleToSite::generatuniqueID()
{
    static int s_id = 0;
    return ++s_id;
}




