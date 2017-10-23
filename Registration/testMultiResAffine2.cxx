#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>

using namespace std;

int main()
{
	double translationFactor[] = {0.05,0.1,0.2,0.3};
	int iterationsNumber[] = {  100 };
	int numberOfLevels[] = {1};
	int numberOfBins[] = { 50 };
	for (int b = 0; b < sizeof(numberOfBins) / sizeof(numberOfBins[0]); b++)
	{
		for (int j = 0; j < sizeof(numberOfLevels) / sizeof(numberOfLevels[0]); j++)
		{
			for (int k = 0; k < sizeof(iterationsNumber) / sizeof(iterationsNumber[0]); k++)
			{
				for (int i = 0; i < sizeof(translationFactor) / sizeof(translationFactor[0]); i++)
				{
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/";
					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20150801/MultiResAffine/Level%d_TransScale%g_IterationNum%d_", numberOfLevels[j],  translationFactor[i], iterationsNumber[k]);

					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);
					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%sdicom_ear_male_leftEar.nii", fixedfolder);


					char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
					char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";

					char outputfolder[300];
					sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/MultiResAffine_first/Level%d_Bin%d_TransScale%g_IterationNum%d_", numberOfLevels[j], numberOfBins[b], translationFactor[i], iterationsNumber[k]);

					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%sfirst_base_line_ear_matched.nii", fixedfolder);
					char movingImageFile[300];
					sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);


					char registeredImageFile[200];
					sprintf_s(registeredImageFile, " %sRegistrationOutput.nii", outputfolder);

					char differenceBeforeFile[200];
					sprintf_s(differenceBeforeFile, " %sdifferenceBeforeRegistration.nii", outputfolder);

					char differenceAfterFile[200];
					sprintf_s(differenceAfterFile, " %sdifferenceAfterRegistration.nii", outputfolder);

					char labelmapOutput_SL[200];
					sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_SM[200];
					sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_ST[200];
					sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_SV[200];
					sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", outputfolder);

					char labelmapInput_SV[200];
					sprintf_s(labelmapInput_SV, " %satlasVolume_SV_labelmap.nii", inputfolder);
					char labelmapInput_ST[200];
					sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", inputfolder);
					char labelmapInput_SM[200];
					sprintf_s(labelmapInput_SM, " %satlasVolume_SM_labelmap.nii", inputfolder);
					char labelmapInput_SL[200];
					sprintf_s(labelmapInput_SL, " %satlasVolume_SL_labelmap.nii", inputfolder);

					//char labelmapInput_SV[300];
					//sprintf_s(labelmapInput_SV, " %sSV_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_ST[300];
					//sprintf_s(labelmapInput_ST, " %sST_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_SM[300];
					//sprintf_s(labelmapInput_SM, " %sSM_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_SL[300];
					//sprintf_s(labelmapInput_SL, " %sSL_labelmap_registrationOutput.nii", inputfolder);

					char parameter[50];
					sprintf_s(parameter, " %d %g %d %d", iterationsNumber[k], translationFactor[i], numberOfLevels[j], numberOfBins[b]);

					char command[3000] = "D:/ITK4.8/ITK_Projects_build/Registration/Debug/MultiResAffine3D2.exe ";

					strcat_s(command, fixedImageFile);
					strcat_s(command, movingImageFile);
					strcat_s(command, registeredImageFile);

					strcat_s(command, differenceBeforeFile);
					strcat_s(command, differenceAfterFile);
					strcat_s(command, parameter);
					/*strcat_s(command, labelmapInput_SV);
					strcat_s(command, labelmapOutput_SV);
					strcat_s(command, labelmapInput_ST);
					strcat_s(command, labelmapOutput_ST);
					strcat_s(command, labelmapInput_SM);
					strcat_s(command, labelmapOutput_SM);
					strcat_s(command, labelmapInput_SL);
					strcat_s(command, labelmapOutput_SL);*/

					cout << command << endl;
					std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

					system(command);
					std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

					std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;
				}
			}
		}
	}
	return 0;
}