#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>
#include <complex>
#include <cassert>
using namespace std;


class QuantumState
{
public:

    QuantumState()
    {
        State = nullptr;
        quantityOfQubits = 0;
    }

    ~QuantumState()
    {
        delete[] State;
    }

    void SetQuantumState (const int quantityOfQubits, const int* State)
    {
//        if (this->State != nullptr)
//        {
//            delete []this->State;
//        }
        this->quantityOfQubits = quantityOfQubits;
        this->State = new int[this->quantityOfQubits];
        for (int j = 0; j < quantityOfQubits; ++j)
        {
            CheckRange(j);
            this->State[j] = State[j];
        }
    }

    void setQuantityOfQubits (const int quantityOfQubits)
    {
        if (this->State != nullptr)
        {
            delete []this->State;
        }
        this->quantityOfQubits = quantityOfQubits;
        this->State = new int[this->quantityOfQubits];

    }

    QuantumState &operator=(const QuantumState&other)
    {
        delete[] this->State;
        this->quantityOfQubits = other.quantityOfQubits;
        this->State = new int[quantityOfQubits];
        for (int j = 0; j < quantityOfQubits; ++j)
        {
            CheckRange(j);
            this->State[j] = other.State[j];
        }
    }

    int &operator[](const int index)
    {
        CheckRange(index);
        return State[index];
    }

    bool operator==(const QuantumState&other)
    {
        if (this->quantityOfQubits != other.quantityOfQubits)
        {
            return false;
        }
        for (int j = 0; j < this->quantityOfQubits; ++j)
        {
            if (this->State[j] != other.State[j])
            {
                return false;
            }
        }
        return true;
    }

    void Sx(const int firstQubit, const int secondQubit)
    {
        CheckRange(firstQubit);
        CheckRange(secondQubit);
        State[firstQubit] = !State[firstQubit];
        State[secondQubit] = !State[secondQubit];
    }

    void Sy(const int firstQubit, const int secondQubit)
    {
        CheckRange(firstQubit);
        CheckRange(secondQubit);
        State[firstQubit] = !State[firstQubit];
        State[secondQubit] = !State[secondQubit];
    }

    void CNOT()
    {
        for (int i = 0; i < (quantityOfQubits - 1); ++i)
        {
            if (State[i])
            {
                State[i + 1] = !State[i + 1];
            }

        }
    }

private:
    int *State;
    int quantityOfQubits;

    void CheckRange(const int index)
    {
        assert(index >=0 && index < quantityOfQubits);
    }

};

class WaveFunction
{
public:

    WaveFunction()
    {
        this->quantityOfStates = 2;
        this->quantityOfQubits = 1;
        this->Function = new QuantumState[this->quantityOfStates];
        this->Coefficients = new complex<double>[this->quantityOfStates];
        int State[1] = {0};
        Function[0].SetQuantumState(quantityOfQubits, State);
        State[0] = 1;
        Function[1].SetQuantumState(quantityOfQubits, State);
        resetCoefficients();
    }

    WaveFunction(const int quantityOfQubits);


    WaveFunction(const int quantityOfQubits, const int stateNumber);

    WaveFunction(const WaveFunction& other)
    {
        this->quantityOfStates = other.quantityOfStates;
        this->quantityOfQubits = other.quantityOfQubits;
        Function = new QuantumState[other.quantityOfStates];
        this->Coefficients = new complex<double>[other.quantityOfStates];
        for (int j = 0; j < quantityOfStates; ++j)
        {
            Function[j] = other.Function[j];
            Coefficients[j] = other.Coefficients[j];
        }
    }

    ~WaveFunction()
    {
        delete[] Function;
        delete[] Coefficients;
    }

    void printFunction()
    {
        for (int j = 0; j < quantityOfStates; ++j)
        {
            CheckRange(j);
            cout << Coefficients[j] << " * |";
            for (int k = (quantityOfQubits - 1); k >= 0; --k)
            {
                cout << " " << Function[j][k];
            }
            cout << " >";
            if (j != quantityOfStates - 1)
            {
                cout << " + ";
            }
        }
        cout << endl;
    }

    void printHamiltonian()
    {
        complex<double> **CC = buildHamiltonian();
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            for (int i = 0; i < this->quantityOfStates; ++i)
            {
                cout << CC[j][i] << " ";
            }
            cout << endl;
        }
    }

    void printHamiltonian(double B, double J)
    {
        complex<double> **CC = buildHamiltonian(B, J);
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            for (int i = 0; i < this->quantityOfStates; ++i)
            {
                cout << CC[j][i] << " ";
            }
            cout << endl;
        }
    }

    int & getQuantityOfStates ()
    {
        return quantityOfStates;
    }

    int & getQuantityOfQubits ()
    {
        return quantityOfQubits;
    }

    complex<double>&operator[] (int index)
    {
        CheckRange(index);
        return Coefficients[index];
    }

    WaveFunction &operator=(const WaveFunction &other)
    {
        delete[] this->Coefficients;
        delete[] this->Function;
        this->quantityOfQubits = other.quantityOfQubits;
        this->quantityOfStates = other.quantityOfStates;
        this->Coefficients = new complex<double>[this->quantityOfStates];
        this->Function = new QuantumState[this->quantityOfStates];
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            this->Coefficients[j] = other.Coefficients[j];
            this->Function[j] = other.Function[j];
        }
        return  *this;
    }

    WaveFunction operator*(const WaveFunction &other)
    {
        WaveFunction temp;
        delete[] temp.Coefficients;
        delete[] temp.Function;
        temp.quantityOfQubits = this->quantityOfQubits + other.quantityOfQubits;
        temp.quantityOfStates = this->quantityOfStates * other.quantityOfStates;
        temp.Coefficients = new complex<double>[temp.quantityOfStates];
        temp.Function = new QuantumState[temp.quantityOfStates];
        for (int i = 0; i < temp.quantityOfStates; ++i)
        {
            temp.Function[i].setQuantityOfQubits(temp.quantityOfQubits);
            temp.Coefficients[i] = 1;
        }
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            for (int i = 0; i < other.quantityOfStates; ++i)
            {
                for (int k = 0; k < this->quantityOfQubits; ++k)
                {
                    temp.Function[i + (j * other.quantityOfStates)][k] = this->Function[j][k];
                }
                for (int l = (this->quantityOfQubits); l < (this->quantityOfQubits + other.quantityOfQubits); ++l)
                {
                    temp.Function[i + (j * other.quantityOfStates)][l] = other.Function[i][l - this->quantityOfQubits];
                }
            }
        }
        return  temp;
    }

    WaveFunction operator+( WaveFunction &other)
    {
        WaveFunction temp(*this);
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            temp.Coefficients[j] = temp.Coefficients[j] + other.Coefficients[j];
        }
        return temp;
    }

    void BSz(const int qubitNumber, const double B, const int j)
    {
        WaveFunction temp(this->quantityOfQubits, j);
        temp.Sz(qubitNumber);
        for (int i = 0; i < this->quantityOfStates; ++i) {
            this->Coefficients[i] += temp.Coefficients[i] * B;
        }
    }

    void S(const int firstQubit, const int secondQubit)
    {
        WaveFunction temp1(*this);
        WaveFunction temp2(*this);
        WaveFunction temp3(*this);
        temp1.Sx(firstQubit, secondQubit);
        temp2.Sy(firstQubit, secondQubit);
        temp3.Sz(firstQubit, secondQubit);
        temp1 = temp1 + temp2;
        temp1 = temp1 + temp3;
        for (int i = 0; i < this->quantityOfStates; ++i) {
            this->Coefficients[i] = temp1.Coefficients[i];
        }
    }

    void CNOT()
    {
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            Function[j].CNOT();
        }
        CheckStates();
    }

    void hamiltonian();

    void hamiltonian(double B, double J);

    void rotationMatrix(const double* qubitAngles);

    void cnot ();

    complex<double> conjugateTranspose(const WaveFunction &other)
    {
        complex<double> temp = {0, 0};
        try
        {
            if (this->quantityOfStates != other.quantityOfStates)
            {
                throw "this->quantityOfStates != other.quantityOfStates";
            }
            for (int i = 0; i < this->quantityOfStates; ++i)
            {
                temp = temp + (this->Coefficients[i] * conj(other.Coefficients[i]));
            }
            return temp;
        }
        catch (char *error)
        {
            cout << error << endl;
        }
    }

private:
    QuantumState *Function;
    complex<double> *Coefficients;
    int quantityOfStates;
    int quantityOfQubits;

    void resetCoefficients()
    {
        for (int j = 0; j < quantityOfStates; ++j)
        {
            CheckRange(j);
            Coefficients[j] = complex<double> (1, 0);
        }
    }

    void resetCoefficients(const int stateNumber)
    {
        for (int j = 0; j < quantityOfStates; ++j)
        {
            CheckRange(j);
            if (stateNumber != -1) {
                CheckRange(stateNumber);

            }
            if (j == stateNumber)
            {
                Coefficients[j] = complex<double> (1, 0);

            }
            else
            {
                Coefficients[j] = complex<double> (0, 0);

            }
        }
    }

    void CheckRange(const int index)
    {
        assert(index >=0 && index < quantityOfStates);
    }

    void CheckStates()
    {
        WaveFunction reference(this->quantityOfQubits, -1);
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            for (int k = 0; k < reference.quantityOfStates; ++k)
            {
                if (this->Function[j] == reference.Function[k])
                {
                    reference.Coefficients[k] += this->Coefficients[j];
                }
            }
        }
        *this = reference;
    }

    void Sx(const int firstQubit, const int secondQubit)
    {
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            Coefficients[j] = 0.25 * Coefficients[j];
            Function[j].Sx(firstQubit, secondQubit);
        }
        CheckStates();
    }

    void Sy(const int firstQubit, const int secondQubit)
    {
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            if (Function[j][firstQubit])
            {
                Coefficients[j] = (complex<double> (0, -0.5)) * Coefficients[j];
            }
            else
            {
                Coefficients[j] = (complex<double> (0, 0.5)) * Coefficients[j];
            }
            if (Function[j][secondQubit])
            {
                Coefficients[j] = (complex<double> (0, -0.5)) * Coefficients[j];
            }
            else
            {
                Coefficients[j] = (complex<double> (0, 0.5)) * Coefficients[j];
            }
            Function[j].Sy(firstQubit, secondQubit);
        }
        CheckStates();
    }

    void Sz(const int firstQubit, const int secondQubit)
    {
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            if (Function[j][firstQubit])
            {
                Coefficients[j] = (complex<double> (0.5, 0)) * Coefficients[j];
            }
            else
            {
                Coefficients[j] = (complex<double> (- 0.5, 0)) * Coefficients[j];
            }
            if (Function[j][secondQubit])
            {
                Coefficients[j] = (complex<double> (0.5, 0)) * Coefficients[j];
            }
            else
            {
                Coefficients[j] = (complex<double> (-0.5, 0)) * Coefficients[j];
            }
        }
        CheckStates();
    }

    void Sz(const int firstQubit)
    {
        for (int j = 0; j < this->quantityOfStates; ++j)
        {
            if (Function[j][firstQubit])
            {
                Coefficients[j] = 0.5 * Coefficients[j];
            }
            else
            {
                Coefficients[j] = -0.5 * Coefficients[j];
            }
        }
        CheckStates();
    }

    complex<double> ** buildHamiltonian ();

    complex<double> ** buildHamiltonian (double B, double J);

    complex<double> **buildRotationMatrix(const double* QubitAngles);

    complex<double> **buildCNOT ();

};

complex<double> ** matrixMuptiplication (complex<double>** A, complex<double>** B, int Arows, int Acolms, int Brows)
{
    complex<double> **CC;
    CC = new complex<double>* [Brows];
    for (int i = 0; i < Brows; ++i)
    {
        CC[i] = new complex<double> [Acolms];
    }
    for (int j = 0; j < Brows; ++j)
    {
        for (int k = 0; k < Acolms; ++k)
        {
            CC[j][k] = 0;
            for (int l = 0; l < Arows; ++l)
            {
                CC[j][k] += B[j][l] * A[l][k];
            }
        }
    }
    for (int k = 0; k < Arows; ++k)
    {
        delete[] A[k];
    }
    delete[]A;
    for (int l = 0; l < Brows; ++l)
    {
        delete[] B[l];

    }
    delete[]B;
    return CC;
}

complex<double> ** matrixMuptiplication (WaveFunction &A, complex<double>** B, int Brows)
{
    complex<double> **CC;
    CC = new complex<double>* [Brows];
    for (int i = 0; i < Brows; ++i)
    {
        CC[i] = new complex<double> [1];
    }
    for (int j = 0; j < A.getQuantityOfStates(); ++j)
    {
        CC[j][0] = 0;
        for (int i = 0; i < Brows; ++i)
        {
            CC[j][0] = CC[j][0] + B[j][i] * A[i];
        }
    }
//    for (int j = 0; j < Brows; ++j)
//    {
//        CC[j][0] = 0;
//        for (int l = 0; l < A.getQuantityOfStates(); ++l)
//        {
//            CC[j][0] += B[j][l] * A[l];
//        }
//    }

    for (int k = 0; k < Brows; ++k)
    {
        delete []B[k];
    }
    delete[]B;
    return CC;
}

complex<double> ** CartesianProduct (complex<double>** A, complex<double>** B, int Arows, int Acolms, int Brows, int Bcolms)
{
    complex<double> **CC = new complex<double>* [Arows  * Brows]() ;
//    CC = new complex<double>* [Arows  * Brows]();
    for (int i = 0; i < Arows  * Brows; ++i)
    {
        CC[i] = new complex<double> [Acolms * Bcolms]();
    }

    for (int j = 0; j < (Arows  * Brows); j = j + Brows)
    {
        for (int k = 0; k < (Acolms * Bcolms); k = k + Bcolms)
        {
            for (int l = 0; l < Brows; ++l)
            {
                for (int m = 0; m < Bcolms; ++m)
                {
                    CC[j + l][k + m] = A[j / Brows][k / Bcolms] * B[l][m];
                }
            }
        }
    }
    return CC;
}

WaveFunction::WaveFunction(const int quantityOfQubits)
{
    this->quantityOfStates = pow(2, quantityOfQubits);
    this->quantityOfQubits = quantityOfQubits;
    this->Function = new QuantumState[this->quantityOfStates];
    this->Coefficients = new complex<double>[this->quantityOfStates];
    WaveFunction temp;
    *this = temp;
    for (int i = 1; i < quantityOfQubits; ++i)
    {
        if (quantityOfQubits != 1)
        {
            *this = temp * *this;
        }
    }
    this->resetCoefficients();
}

WaveFunction::WaveFunction(const int quantityOfQubits, const int stateNumber)
{
    this->quantityOfStates = pow(2, quantityOfQubits);
    this->quantityOfQubits = quantityOfQubits;
    this->Function = new QuantumState[this->quantityOfStates];
    this->Coefficients = new complex<double>[this->quantityOfStates];
    WaveFunction temp;
    *this = temp;
    for (int i = 1; i < quantityOfQubits; ++i)
    {
        if (quantityOfQubits != 1)
        {
            *this = temp * *this;
        }
    }
    this->resetCoefficients(stateNumber);
}

complex<double> ** WaveFunction::buildHamiltonian ()
{
    complex<double> **CC;
    CC = new complex<double>* [this->getQuantityOfStates()];
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        CC[i] = new complex<double> [this->getQuantityOfStates()];
    }
    switch (this->getQuantityOfStates())
    {
        case 2:
        {
            for (int j = 0; j < this->getQuantityOfStates(); ++j)
            {
                WaveFunction temp(this->getQuantityOfQubits(), j);
                temp.Sz(0);
                for (int k = 0; k < this->getQuantityOfStates(); ++k)
                {
                    CC[j][k] = temp[k];
                }

            }
            break;
        }
        default:
        {
            for (int j = 0; j < this->getQuantityOfStates(); ++j)
            {
                WaveFunction temp(this->getQuantityOfQubits(), j);

                for (int i = 1; i < this->quantityOfQubits; ++i) {
                    temp.S(i - 1 , i);

                }
                if (this->getQuantityOfQubits() != 2)
                {
                    temp.S(this->getQuantityOfQubits() - 1, 0);
                }

                for (int k = 0; k < this->getQuantityOfStates(); ++k)
                {
                    CC[j][k] = temp[k];
                }

            }
        }
    }
    return CC;
}

complex<double> ** WaveFunction::buildHamiltonian (double B, double J)
{
    complex<double> **CC;
    CC = new complex<double>* [this->getQuantityOfStates()];
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        CC[i] = new complex<double> [this->getQuantityOfStates()];
    }
    switch (this->getQuantityOfStates())
    {
        case 2:
        {
            for (int j = 0; j < this->getQuantityOfStates(); ++j)
            {
                WaveFunction temp(this->getQuantityOfQubits(), j);
                temp.Sz(0);
                for (int k = 0; k < this->getQuantityOfStates(); ++k)
                {
                    CC[j][k] = temp[k];
                }

            }
            break;
        }
        default:
        {
            for (int j = 0; j < this->getQuantityOfStates(); ++j)
            {
                WaveFunction temp(this->getQuantityOfQubits(), j);
                for (int i = 1; i < this->quantityOfQubits; ++i) {
                    temp.S(i - 1 , i);
                }
                if (this->getQuantityOfQubits() != 2)
                {
                    temp.S(this->getQuantityOfQubits() - 1, 0);
                }
                for (int l = 0; l < this->quantityOfQubits; ++l) {
                    temp.BSz(l, B, j);
                }
                for (int k = 0; k < this->getQuantityOfStates(); ++k)
                {
                    CC[j][k] = J * temp[k];
                }

            }
        }
    }
    return CC;
}

void WaveFunction::hamiltonian()
{
    complex<double>** temp = nullptr;
    temp = matrixMuptiplication(*this, this->buildHamiltonian(), this->getQuantityOfStates());
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        this->Coefficients[i]  = temp[i][0];
    }
    for (int j = 0; j < this->getQuantityOfStates(); ++j) {
        delete []temp[j];
    }
    delete []temp;
}

void WaveFunction::hamiltonian(double B, double J)
{
    complex<double>** temp = nullptr;
    temp = matrixMuptiplication(*this, this->buildHamiltonian(B, J), this->getQuantityOfStates());
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        this->Coefficients[i]  = temp[i][0];
    }
    for (int j = 0; j < this->getQuantityOfStates(); ++j) {
        delete []temp[j];
    }
    delete []temp;
}

complex<double>** WaveFunction::buildRotationMatrix(const double* QubitAngles)
{
    complex<double> **CC;
    CC = new complex<double>* [this->getQuantityOfStates()];
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        CC[i] = new complex<double> [this->getQuantityOfStates()];
    }

    for (int l = 0; l < this->getQuantityOfStates(); ++l) {
        for (int i = 0; i < this->getQuantityOfStates(); ++i) {
            CC[l][i] = 0;
        }
    }

//    cout << "Angles" << endl;
//    cout << QubitAngles  << " "<< QubitAngles[1] << " "<< QubitAngles[2];
//    cout << endl;

    CC[0][0] = complex<double>(cos(QubitAngles[1] / 2), 0);
    CC[0][1] = complex<double>( - sin(QubitAngles[1] / 2) * cos(QubitAngles[2]), - sin(QubitAngles[1] / 2) * sin(QubitAngles[2]));
    CC[1][0] = complex<double>( sin(QubitAngles[1] / 2) * cos(QubitAngles[0] ), sin(QubitAngles[1] / 2) * sin(QubitAngles[0]  ));
    CC[1][1] = complex<double>( cos(QubitAngles[1] / 2) * cos( QubitAngles[0] + QubitAngles[2]), cos(QubitAngles[1] / 2) * sin( QubitAngles[0] + QubitAngles[2]));

    int CCrows = 2;
    int CCcolms = 2;

//    for (int j = 0; j < CCrows; ++j) {
//        for (int i = 0; i < CCcolms; ++i) {
//            cout << CC[j][i] << " ";
//        }
//        cout << endl;
//    }

    for (int i = 1; i < this->getQuantityOfQubits(); ++i)
    {
        if (this->getQuantityOfQubits() != 1)
        {
//            cout << "Angles" << endl;
//            cout << QubitAngles[(i * 3)] << " "<< QubitAngles[1 + (i * 3)] << " "<< QubitAngles[2 + (i * 3)];
//            cout << endl;

            complex<double> **C1;
            C1 = new complex<double>* [2];
            for (int i = 0; i < 2; ++i)
            {
                C1[i] = new complex<double> [2];
            }

            C1[0][0] = complex<double>(cos(QubitAngles[1 + (i * 3)] / 2), 0);
            C1[0][1] = complex<double>( - sin(QubitAngles[1 + (i * 3)] / 2) * cos(QubitAngles[2 + (i * 3)]), - sin(QubitAngles[1 + (i * 3)] / 2) * sin(QubitAngles[2 + (i * 3)]));
            C1[1][0] = complex<double>( sin(QubitAngles[1 + (i * 3)] / 2) * cos(QubitAngles[0 + (i * 3)]), sin(QubitAngles[1 + (i * 3)] / 2) * sin(QubitAngles[0 + (i * 3)]));
            C1[1][1] = complex<double>( cos(QubitAngles[1 + (i * 3)] / 2) * cos( QubitAngles [0 + (i * 3)] + QubitAngles[2 + (i * 3)]), cos(QubitAngles[1 + (i * 3)] / 2) * sin( QubitAngles [0 + (i * 3)] + QubitAngles[2 + (i * 3)]));
//            for (int j = 0; j < 2; ++j) {
//                for (int i = 0; i < 2; ++i) {
//                    cout << C1[j][i] << " ";
//                }
//                cout << endl;
//            }
//            cout << "++++++++++++++++++++++++++++" << endl;

//            CC = CartesianProduct(C1, CC, 2, 2, CCrows, CCcolms);
            CC = CartesianProduct(CC, C1, CCrows, CCcolms, 2, 2);
            CCrows *= 2;
            CCcolms *= 2;

//            for (int j = 0; j < CCrows; ++j) {
//                for (int i = 0; i < CCcolms; ++i) {
//                    cout << CC[j][i] << " ";
//                }
//                cout << endl;
//            }
//            cout << "=================" << endl;

            for (int k = 0; k < 2; ++k)
            {
                delete []C1[k];
            }
            delete[]C1;
        }
    }
//    for (int j = 0; j < CCrows ; ++j) {
//        for (int i = 0; i < CCcolms; ++i) {
//            cout << CC[j][i] << " ";
//        }
//        cout << endl;
//    }




    return CC;
}

void WaveFunction::rotationMatrix(const double *qubitAngles)
{
    complex<double>** temp;
    temp = matrixMuptiplication(*this, buildRotationMatrix(qubitAngles), this->getQuantityOfStates());
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        this->Coefficients[i]  = temp[i][0];
    }
    for (int j = 0; j < this->getQuantityOfStates(); ++j) {
        delete []temp[j];
    }
    delete []temp;
}

complex<double> ** WaveFunction::buildCNOT()
{
    assert(this->quantityOfQubits > 1);

    complex<double> **CC;
    CC = new complex<double>* [this->getQuantityOfStates()];
    for (int i = 0; i < this->getQuantityOfStates(); ++i)
    {
        CC[i] = new complex<double> [this->getQuantityOfStates()];
    }

    for (int j = 0; j < this->getQuantityOfStates(); ++j)
    {
        WaveFunction temp(this->getQuantityOfQubits(), j);
        temp.CNOT();
        for (int k = 0; k < this->getQuantityOfStates(); ++k)
        {
            CC[k][j] = temp[k];
        }
    }
    return CC;
}

void WaveFunction::cnot()
{
    if (this->quantityOfQubits != 1)
    {
        complex<double>** temp = nullptr;
        temp = matrixMuptiplication(*this, buildCNOT(), this->getQuantityOfStates());
        for (int i = 0; i < this->getQuantityOfStates(); ++i)
        {
            this->Coefficients[i]  = temp[i][0];
        }
        for (int j = 0; j < this->getQuantityOfStates(); ++j) {
            delete []temp[j];
        }
        delete []temp;
    }
}

void ResetValues( double *Angles, int numberOfAngles)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist1(0, 200 * M_PI);
    uniform_int_distribution<> dist2(0, 100 * M_PI);

    for (int j = 0; j < numberOfAngles; j = j + 3)
    {
        Angles [j] = (double) dist1(gen) / 100;
    }
    for (int k = 1; k < numberOfAngles; k = k + 3)
    {
        Angles [k] = (double) dist2(gen) / 100;

    }
    for (int l = 2; l < numberOfAngles; l = l + 3)
    {
        Angles [l] = (double) dist1(gen) / 100;

    }
}

void  NewBernDistrib(int*  BernDistrib, int numberOfAngles)
{
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0,1);
    for (int j = 0; j < numberOfAngles; ++j)
    {
//        BernDistrib [j] = 1;
        BernDistrib [j] = dist(gen);
    }
}

double Energy (const double *qubitAngles, int quantityOfQubits, int d, double B, double J)
{
    double EnergyValue;
    WaveFunction Psi(quantityOfQubits, 0);
    for (int j = 0; j < d; ++j)
    {
        double *QubitAnglesD  = new double [quantityOfQubits * 3];
        for (int i = 0; i < quantityOfQubits * 3; ++i) {
            QubitAnglesD[i] = qubitAngles[i + (j * quantityOfQubits * 3)];
        }
        Psi.rotationMatrix(QubitAnglesD);
        Psi.CNOT();
        delete [] QubitAnglesD;
    }
    WaveFunction PsiCong(Psi);
    Psi.hamiltonian(B, J);
//    cout << Psi.conjugateTranspose(PsiCong) << endl;
    EnergyValue = Psi.conjugateTranspose(PsiCong).real();
    return EnergyValue;
}

void AngleRestriction(double *Angles, int quantityOfAngles)
{
    for (int k = 1; k < quantityOfAngles; k = k + 3)
    {
        if (Angles[k] > M_PI || Angles[k] < 0) {
            bool Restriction = true;
            while (Restriction) {
                if (Angles[k] >= 2 * M_PI)
                {
                    Angles[k] = Angles[k] - 2 * M_PI;
                }
                if (Angles[k] > M_PI && Angles[k]< 2 * M_PI) {
                    Angles[k] = 2 * M_PI - Angles[k];
                    Angles[k + 1] = Angles[k - 1] + M_PI;
                }
                if (Angles[k] <= -2 * M_PI) {
                    Angles[k] = Angles[k] + 2 * M_PI;
                }
                if (Angles[k] > -2 * M_PI && Angles[k]< 0) {
                    Angles[k] = -Angles[k];
                    Angles[k + 1] = Angles[k - 1] + M_PI;
                }
                if (Angles[k] >= 0 && Angles[k] <= M_PI) {
                    Restriction = false;
                }
            }
        }
    }
    for (int j = 0; j < quantityOfAngles; j = j + 3)
    {
        if (Angles[j] > 2 * M_PI || Angles  [j] < 0) {
            bool Restrition = true;
            while (Restrition) {
                if (Angles[j] > 2 * M_PI) {
                    Angles[j] = Angles  [j] - 2 * M_PI;
                }
                if (Angles[j]< 0) {
                    Angles[j] = Angles  [j] + 2 * M_PI;
                }
                if (Angles[j] >= 0 && Angles[j] <= 2 * M_PI) {
                    Restrition = false;
                }
            }
        }
    }
    for (int l = 2; l < quantityOfAngles; l = l + 3)
    {
        if (Angles[l] > 2 * M_PI || Angles  [l] < 0) {
            bool Restrition = true;
            while (Restrition) {
                if (Angles[l] > 2 * M_PI) {
                    Angles[l] = Angles  [l] - 2 * M_PI;
                }
                if (Angles[l]< 0) {
                    Angles[l] = Angles  [l] + 2 * M_PI;
                }
                if (Angles[l] >= 0 && Angles[l] <= 2 * M_PI) {
                    Restrition = false;
                }
            }
        }
    }
}

int main()
{
    int quantityOfQubits = 2;
    int d = 1;
    int ll = 0;
    int kk = 0;
    double C = 0.1;
    double Alpha = 0.602;
    double  Gamma = 0.101;
    double J = 1;
    int number_of_iteration = 500;
    double B = 0;
    int quantityOfAngles = quantityOfQubits * 3 * d;

//    double *Angles = new double[quantityOfAngles]{3.14159, 2.4, 5.90593, 5.44565, 1.02, 3.70057 };
////    double *Angles = new double[quantityOfAngles]{1.47278, 2.31286, 5.15356, 3.14176, 1.5708, 2.42552, 5.8797, 3.1414, 0.540039 };
////    ResetValues(Angles, quantityOfAngles);
//    for (int i = 0; i < quantityOfAngles; i += 3)
//    {
//        Angles[i] = 0;
////        Angles[i + 1] = 0;
//        Angles[i + 2] = 0;
//    }
//
//
//    for (int j = 0; j < quantityOfAngles; ++j) {
//        cout << Angles[j] << " ";
//    }
//    cout << endl;
//    WaveFunction Function(quantityOfQubits, 0);
//    Function.printFunction();
//    Function.rotationMatrix(Angles);
////    Function.CNOT();
//    Function.printFunction();
//    cout << Energy(Angles, quantityOfQubits, 1, B);
//    delete[] Angles;

//
    double *Angles = new double[quantityOfAngles];
    double *AnglesTemp = new double[quantityOfAngles];
    double *AnglesPlus = new double[quantityOfAngles];
    double *AnglesMinus = new double[quantityOfAngles];
    int *BernDistrib = new int[quantityOfAngles];
    ofstream Energy_Iteration("Energy_Iteration.txt", ios_base::out | ios_base::trunc);
    ofstream Energy_MagneticField("Energy_MagneticField.txt", ios_base::out | ios_base::trunc);
//    for (B = 0; B < 2; B = B + 0.1) {
    ResetValues(Angles, quantityOfAngles);
    double E = Energy(Angles, quantityOfQubits, d, B, J);
    for (int i = 1; i < number_of_iteration; ++i) {
        Energy_Iteration << i << " " << E  << endl;
//            cout << i << " " << E << endl;
        NewBernDistrib(BernDistrib, quantityOfAngles);
        double Ci = C / pow(i, Gamma);
        for (int j = 0; j < quantityOfAngles; ++j)
        {
            AnglesPlus[j] = Angles[j] + (Ci * BernDistrib[j]);
            AnglesMinus[j] = Angles[j] - (Ci * BernDistrib[j]);
        }
        double *Gradient = new double [quantityOfAngles];
        double gradient_k = ((Energy(AnglesPlus, quantityOfQubits, d, B, J) - Energy(AnglesMinus, quantityOfQubits, d, B, J)) / (2 * Ci));
        for (int k = 0; k < quantityOfAngles; ++k)
        {
            Gradient[k] = gradient_k * BernDistrib[k];
        }
        double A = C;
//            if (i == 1)
//            {
//                double *AnglesBuffer = new double [NumberOfQubits * 3];
//                double Energy1 = 0;
//                AnglesBuffer = Angles;
//                for (int k = 0; k < 25; ++k) {
//                    NewBernDistrib();
//                    for (int j = 0; j < Angles[1][j]; ++j) {
//                        Angles[1][j] = AnglesBuffer[j] + Ci * BernDistrib[j];
//                        Angles[2][j] = AnglesBuffer[j] - Ci * BernDistrib[j];
//                    }
//                    Energy1 = Energy1 + (Energy(Angles[1], quantityOfQubits, d) - Energy(Angles[2], quantityOfQubits, d));
//                }
//                delete[] AnglesBuffer;
//                A = (2 * M_PI / 5) * C / (Energy1 / 25);
//            }
        double Ai = A / pow(i, Alpha);
        Ai = C;
        for (int l = 0; l < quantityOfAngles; ++l) {
            AnglesTemp[l] = Angles[l] - Ai * Gradient[l];
        }
        delete[] Gradient;
        double Etemp = Energy(AnglesTemp, quantityOfQubits, d, B, J);
        if (Etemp > E) {
//                continue;
        }
        E = Etemp;
        Angles = AnglesTemp;
//            AngleRestriction(Angles, quantityOfAngles);
    }
    cout << E << endl;
//        if (E < -2.8)
//        {
//            ll++;
//        }
//        if (E > -2.8)
//        {
//            kk++;
//        }
//    for (int m = 0; m < quantityOfAngles; ++m) {
//        cout << Angles[m] << " ";
//    }
    Energy_MagneticField << B << "\t " << E << "\t ";
    for (int m = 0; m < quantityOfAngles; ++m)
    {
        Energy_MagneticField << Angles[m] << " ";
    }
    Energy_MagneticField << endl;
//    }
//    cout << "ll: " << ll << " kk: " << kk;
    for (int n = 0; n < quantityOfAngles; ++n) {
        cout << Angles[n] << " ";
    }

    delete[] Angles;
//    delete[] AnglesTemp;
    delete[] AnglesPlus;
    delete[] AnglesMinus;
    delete[] BernDistrib;
    Energy_Iteration.close();
    Energy_MagneticField.close();
    return 0;
}