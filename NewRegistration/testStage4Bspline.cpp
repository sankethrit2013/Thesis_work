#include <iostream>
#include <string>
#include <chrono>
#include "stdio.h"
using namespace std;

int main(int argc, char *argv[])
{
	int numberOfGridNodesInOneDimension[] = { 10 };
	int MaximumNumberOfIterations[] = { 10 };
	for (size_t t = 0; t < 5; t++)
	{
		if (t == 1)
		{
			for (int j = 0; j < sizeof(MaximumNumberOfIterations) / sizeof(MaximumNumberOfIterations[0]); j++)
			{
				for (int k = 0; k < sizeof(numberOfGridNodesInOneDimension) / sizeof(numberOfGridNodesInOneDimension[0]); k++)
				{
					char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/";
					char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/";

					char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/Bspline/";

					char labelmapOutputFolder[300];
					sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/FinalBspline/labelmap_BsplineStage4_Grid%d_Iter%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					char outputfolder[300];
					sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/Stage4_BsplineRegistration_Grid%d_Iterations%d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					/*-----------------------------------------------------------------------------------------------------*/


					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%sscan1.nii", fixedfolder);
					char movingImageFile[300];
					sprintf_s(movingImageFile, "%sBSpline_GridNodeinOneDimension10_Iterations100RegistrationOutput.nii", fixedfolder);

					char maskFile[300];
					sprintf_s(maskFile, "%sST_mask.nii", inputfolder);

					char registeredImageFile[300];
					sprintf_s(registeredImageFile, "%sRegistrationOutput.nii", outputfolder);

					char differenceBeforeFile[300];
					sprintf_s(differenceBeforeFile, "%sdifferenceBeforeRegistration.nii.gz", outputfolder);

					char differenceAfterFile[300];
					sprintf_s(differenceAfterFile, "%sdifferenceAfterRegistration.nii.gz", outputfolder);

					char deformationField[300];
					sprintf_s(deformationField, "%sdeformationField.nrrd", outputfolder);
					/*-----------------------------------------------------------------------------------------------------*/

#pragma region labelMapFileName
					//output
					char labelmapOutput_SL[200];
					sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", labelmapOutputFolder);

					char labelmapOutput_SM[200];
					sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", labelmapOutputFolder);

					char labelmapOutput_ST[200];
					sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", labelmapOutputFolder);

					char labelmapOutput_SV[200];
					sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", labelmapOutputFolder);

					char labelmapOutput_CA[300];
					sprintf_s(labelmapOutput_CA, " %sCA_labelmap_registrationOutput.nii", labelmapOutputFolder);

					//Input
					//Sanketh - change input folder to labelmapInputFolder
					char labelmapInput_SV[200];
					sprintf_s(labelmapInput_SV, " %slabelmap_Deformable_Grid10_Iter100SV_labelmap_registrationOutput.nii", labelmapInputFolder);

					char labelmapInput_ST[200];
					sprintf_s(labelmapInput_ST, " %slabelmap_Deformable_Grid10_Iter100ST_labelmap_registrationOutput.nii", labelmapInputFolder);

					char labelmapInput_SM[200];
					sprintf_s(labelmapInput_SM, " %slabelmap_Deformable_Grid10_Iter100SM_labelmap_registrationOutput.nii", labelmapInputFolder);

					char labelmapInput_SL[200];
					sprintf_s(labelmapInput_SL, " %slabelmap_Deformable_Grid10_Iter100SL_labelmap_registrationOutput.nii", labelmapInputFolder);

					char labelmapInput_CA[300];
					sprintf_s(labelmapInput_CA, " %slabelmap_Deformable_Grid10_Iter100CA_labelmap_registrationOutput.nii", labelmapInputFolder);

#pragma endregion

					char parameter[50];
					sprintf_s(parameter, " %d %d", numberOfGridNodesInOneDimension[k], MaximumNumberOfIterations[j]);

					char command[8000] = "";
					strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/NewRegistration/Debug/Stage4Bspline.exe ");

					strcat_s(command, fixedImageFile);
					strcat_s(command, movingImageFile);
					strcat_s(command, registeredImageFile);
					strcat_s(command, differenceAfterFile);
					strcat_s(command, differenceBeforeFile);
					strcat_s(command, maskFile);
					strcat_s(command, parameter);
					strcat_s(command, deformationField);
					strcat_s(command, labelmapInput_SV);
					strcat_s(command, labelmapOutput_SV);
					strcat_s(command, labelmapInput_ST);
					strcat_s(command, labelmapOutput_ST);
					strcat_s(command, labelmapInput_SM);
					strcat_s(command, labelmapOutput_SM);
					strcat_s(command, labelmapInput_SL);
					strcat_s(command, labelmapOutput_SL);
					strcat_s(command, labelmapInput_CA);
					strcat_s(command, labelmapOutput_CA);

					cout << command << endl;
					std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

					system(command);
					std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

					std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;

				}
			}

		}
	}
	return EXIT_SUCCESS;
}