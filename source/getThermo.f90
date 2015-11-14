! getThermo - A program to calculate Partition Functions and Gibbs Free Energy from selected Normal Modes from GAMESS
! Grupo de Espectroscopia Teórica e Modelagem Molecular - GETMM (DFq/IQ/UFRJ)

FUNCTION QVib(VibModeList, NVibMode, Temperature)
! QVib calculate the Vibrational Partition Function from selected Normal Modes
! Variable Declarations
	REAL(KIND=8), INTENT(IN), DIMENSION(NVibMode) :: VibModeList
	REAL(KIND=8), DIMENSION(NVibMode) :: QVibList
	REAL(KIND=8), PARAMETER :: CPlanck = 6.6260695729E-34						! Planck Constant (J.s)
	REAL(KIND=8), PARAMETER :: LightSpeed = 2.99792458E10						! Speed of Light (cm/s)
	REAL(KIND=8), PARAMETER :: CBoltzmann = 1.380648813E-23						! Boltzmann Constant (J/K)
	REAL(KIND=8) :: QVib, Constant, Temperature
	INTEGER :: NVibMode, Num

	Constant=(CPlanck*LightSpeed)/(CBoltzmann*Temperature)						! Calculate constant part of equation

! Loop to calculate the Vibrational Partition Function to each Normal Mode and the total Vibrational Partition Function
	QVib=1.0
	DO Num = 1, NVibMode
		QVibList(Num)=1.0/(1.0-EXP(-Constant*VibModeList(Num)))					! Calculate individual Partition Functions
		QVib=QVib*QVibList(Num)													! Calculate total Vibrational Partition Function
	END DO
END FUNCTION QVib



FUNCTION QTrans(AMass, NAtoms, Temperature)
! QTrans calculate the Translational Partition Function from selected Temperature
! Variable Declarations
	REAL(KIND=8), INTENT(IN), DIMENSION(NAtoms) :: AMass
	REAL(KIND=8), PARAMETER :: CPlanck = 6.6260695729E-34						! Planck Constant (J.s)
	REAL(KIND=8), PARAMETER :: CBoltzmann = 1.380648813E-23						! Boltzmann Constant (J/K)
	REAL(KIND=8), PARAMETER :: Pressure = 1.01325E5								! Pressure (Pa)
	REAL(KIND=8), PARAMETER :: NAvogadro = 6.0221412927E23						! Avogadro's Number (1/mol)
	REAL(KIND=8), PARAMETER :: Pi = 3.141592653									! Pi number
	REAL(KIND=8) :: Constant
	REAL(KIND=8) :: QTrans, Temperature
	INTEGER :: NAtoms, Num

	Constant=((2.0*Pi/(CPlanck**2.0))**(3.0/2.0))*(CBoltzmann**(5.0/2.0))		! Calculate constant part of equation

! Loop to calculate the total Mass Weight
	TMass=0.0
	DO Num = 1, NAtoms
		TMass=TMass+AMass(Num)													! Calculate individual Mass Weight
	END DO

	TMass=TMass/1000.0															! Total mass weight in International System of Units

	QTrans=Constant*((TMass/NAvogadro)**(3.0/2.0))*(Temperature**(5.0/2.0))/Pressure
END FUNCTION QTrans



FUNCTION QRot(Inertia, Temperature)
! QTrans calculate the Rotational Partition Function from selected Temperature
! Variable Declarations
	REAL(KIND=8), INTENT(IN), DIMENSION(3) :: Inertia
	REAL(KIND=8), PARAMETER :: CPlanck = 6.6260695729E-34						! Planck Constant (J.s)
	REAL(KIND=8), PARAMETER :: CBoltzmann = 1.380648813E-23						! Boltzmann Constant (J/K)
	REAL(KIND=8), PARAMETER :: NAvogadro = 6.0221412927E23						! Avogadro's Number (1/mol)
	REAL(KIND=8), PARAMETER :: Bohr = 5.29177208E-11							! Unit factor conversion from Bohr to Meters
	REAL(KIND=8), PARAMETER :: AMU = 1.0E-3/NAvogadro							! Unit factor conversion from AMU to Kilogram
	REAL(KIND=8), PARAMETER :: Pi = 3.141592653									! Pi number
	REAL(KIND=8), PARAMETER :: NSymmetry = 1.0									! Symmetry Number
	REAL(KIND=8), DIMENSION(3) :: Teta
	REAL(KIND=8) :: QRot, Temperature
	INTEGER :: Num

! Calculate characteristic rotational temperature
	DO Num = 1, 3
		Teta(Num)=(CPlanck**2.0)/(8.0*(Pi**2.0)*Inertia(Num)*AMU*(Bohr**2.0)*CBoltzmann)
	END DO

	QRot=(SQRT(Pi)*(Temperature**(3.0/2.0)))/(NSymmetry*SQRT(Teta(1))*SQRT(Teta(2))*SQRT(Teta(3)))
END FUNCTION QRot



FUNCTION ZPE(VibModeList, NVibMode)
! ZPE calculate the Zero-Point Energy to correct the energy obtained by Harmonic Oscillator Aproximation
! Variable Declarations
	REAL(KIND=8), INTENT(IN), DIMENSION(NVibMode) :: VibModeList
	REAL(KIND=8), DIMENSION(NVibMode) :: ZPEList
	REAL(KIND=8), PARAMETER :: CPlanck = 6.6260695729E-34						! Planck Constant (J.s)
	REAL(KIND=8), PARAMETER :: LightSpeed = 2.99792458E10						! Speed of Light (cm/s)
	REAL(KIND=8), PARAMETER :: NAvogadro = 6.0221412927E23						! Avogadro's Number (1/mol)
	INTEGER :: NVibMode, Num
	REAL(KIND=8) :: ZPE, Constant

	Constant=(CPlanck*LightSpeed*NAvogadro)										! Calculate constant part of equation

! Loop to calculate the Zero-Point Energy to each Normal Mode and the total Zero-Point Energy
	ZPE=0.0
	DO Num = 1, NVibMode
		ZPEList(Num)=(Constant*VibModeList(Num)/2.0)							! Calculate individual Zero-Point Energy
		ZPE=ZPE+ZPEList(Num)													! Calculate total Zero-Point Energy
	END DO
END FUNCTION ZPE



FUNCTION GEnergy(Q, Temperature)
! GEnergy calculate the Gibbs Free Energy contribution to a kind of mode
! Variable Declarations
	REAL(KIND=8), INTENT(IN) :: Q
	REAL(KIND=8), PARAMETER :: CGas = 8.314462175								! Gas Constant (J/K.mol)
	REAL(KIND=8) :: GEnergy, Constant, Temperature

	Constant=(CGas*Temperature)													! Calculate constant part of equation
	GEnergy=-Constant*LOG(Q)													! Calculate and return the Gibbs Free Energy
END FUNCTION GEnergy



FUNCTION GEnergyCorr(QElec, QTrans, QRot, QVib, ZPE, Temperature)
! GEnergyCorr calculate the total Thermal Gibbs Free Energy
! Variable Declarations
	REAL(KIND=8), INTENT(IN) :: QElec, QTrans, QRot, QVib, ZPE
	REAL(KIND=8) :: GEnergy, GEnergyCorr, Temperature

! Calculate and returns Thermal Gibbs Free Energy using Free Energies from Eletronic, Translational, Rotational and Vibrational 
! contributuions with Zero-Point Energy correction
	GEnergyCorr=GEnergy(QElec, Temperature)+GEnergy(QTrans, Temperature)+GEnergy(QRot, Temperature)+GEnergy(QVib, Temperature)+ZPE
END FUNCTION GEnergyCorr



FUNCTION TGEnergy(EEnergy, QElec, QTrans, QRot, QVib, ZPE, Temperature)
! TGEnergy calculate the total Gibbs Free Energy
! Variable Declarations
	REAL(KIND=8), INTENT(IN) :: EEnergy, QElec, QTrans, QRot, QVib, ZPE
	REAL(KIND=8), PARAMETER :: Joule = 4.3597439422E-18							! Energy Unit (J/Hartree)
	REAL(KIND=8), PARAMETER :: NAvogadro = 6.0221412927E23						! Avogadro's Number (1/mol)
	REAL(KIND=8) :: GEnergy, TGEnergy, Temperature

! Calculate and returns Gibbs Free Energy using Free Energies from Eletronic, Translational, Rotational and Vibrational 
! contributuions with Zero-Point Energy correction and Eletronic Energy
	TGEnergy=GEnergy(QElec, Temperature)+GEnergy(QTrans, Temperature)+GEnergy(QRot, Temperature)+GEnergy(QVib, Temperature)+ZPE+(EEnergy*Joule*NAvogadro)
END FUNCTION TGEnergy


PROGRAM getThermo
! Main Program getThermo
! Variable Declarations
	IMPLICIT NONE
	INTEGER :: TMode, NVibMode, Num, NTemp, NAtoms
	INTEGER, ALLOCATABLE, DIMENSION(:) :: IMode
	REAL(KIND=8), DIMENSION(3) :: Inertia
	REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: TModeList, VibModeList, TempList, AMass
	REAL(KIND=8) :: EEnergy, QElec, QTrans, QRot, QVib, ZPE, GEnergy, GEnergyCorr, TGEnergy, TMass
	CHARACTER(100) :: Job, Input

	CALL GETARG(1,Job)															! Get the external argument contents job name
	CALL GETARG(2,Input)														! Get the external argument contents input filename
	OPEN (1, FILE=Input, STATUS='OLD')											! Open input file
	READ (1, *) NAtoms, TMode, NVibMode											! Get number of atoms and of total and vibrational frequencies 

	ALLOCATE(AMass(NAtoms))														! Defines array lenght to contains mass weight from each atom
	ALLOCATE(TModeList(TMode))													! Defines array lenght to contains frequencies
	ALLOCATE(IMode(NVibMode))													! Defines array lenght to contains selected frequencies index
	ALLOCATE(VibModeList(NVibMode))												! Defines array lenght to contais selected frequencies

	READ (1, *) IMode, EEnergy, QElec, AMass, Inertia, TModeList				! Get remaining variables

	READ (1, *) NTemp															! Get number of total temperatures 
	ALLOCATE(TempList(NTemp))													! Defines array lenght to contais selected temperatures
	READ (1, *) TempList														! Get selected temperatures

! Loop to tranfers selected frequencies to a specific array
	DO Num = 1, NVibMode
		VibModeList(Num)=TModeList(Num+TMode-NVibMode)
	END DO

! Header of getThermo Output 
	PRINT *, "###########################################################"
	PRINT *, "#                        getThermo                        #"
	PRINT *, "#         Thermodynamics with selected normal modes       #"
	PRINT *, "###########################################################"

! Print calculation details from input file
	PRINT '(/1X, "Job Name:" T12, A)', Job
	PRINT '(1X, "Total number of normal modes:" T38, I5)', TMode
	PRINT '(1X, "Number of selected normal modes:" T38, I5)', NVibMode
	PRINT '(1X, "Number of selected temperatures:" T38, I5)', NTemp

! Calculate and print Themodynamics Functions
	DO Num=1, NTemp
    	PRINT '(/1X, "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #")'
		PRINT '(1X, "Temperature:" T41, F12.2, 1X, "K")', TempList(Num)
		PRINT '(1X, "Vibrational Partition Function:" T38, Es15.8)', QVib(VibModeList, NVibMode, TempList(Num))
		PRINT '(1X, "Rotational Partition Function:" T38, Es15.8)', QRot(Inertia, TempList(Num))
		PRINT '(1X, "Translational Partition Function:" T38, Es15.8)', QTrans(AMass, NAtoms, TempList(Num))

		PRINT '(/1X, "Vibrational Gibbs Free Energy: " T38, F15.8, 1X, "kJ/mol")', (GEnergy(QVib(VibModeList, NVibMode, TempList(Num)), TempList(Num))+ZPE(VibModeList, NVibMode))/1000.0
		PRINT '(1X, "Rotational Gibbs Free Energy: " T38, F15.8, 1X, "kJ/mol")', (GEnergy(QRot(Inertia, TempList(Num)), TempList(Num)))/1000.0
		PRINT '(1X, "Translational Gibbs Free Energy: " T38, F15.8, 1X, "kJ/mol")', (GEnergy(QTrans(AMass, NAtoms, TempList(Num)), TempList(Num)))/1000.0

        PRINT '(/1X, "Total Thermal Correction:", T38, F15.8, 1X, "kJ/mol")', (GEnergyCorr(QElec, QTrans(AMass, NAtoms, TempList(Num)), QRot(Inertia, TempList(Num)), QVib(VibModeList, NVibMode, TempList(Num)), ZPE(VibModeList, NVibMode), TempList(Num)))/1000.0
		PRINT '(1X, "Total Gibbs Free Energy:", T38, Es15.8, 1X, "kJ/mol")', (TGEnergy(EEnergy, QElec, QTrans(AMass, NAtoms, TempList(Num)), QRot(Inertia, TempList(Num)), QVib(VibModeList, NVibMode, TempList(Num)), ZPE(VibModeList, NVibMode), TempList(Num)))/1000.0
	END DO

END PROGRAM getThermo
