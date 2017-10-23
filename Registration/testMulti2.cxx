#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>

using namespace std;

int main()
{
	double STranslationFactor[] = {2};
	double SRescaleFactor[] = {5}; // changed from 5 to 3
	int iterationsNumber1[] = {200}; //Original 200
	int iterationsNumber2[] = { 100 }; //Original 100
	int iterationsNumber3 = 100 ; //Original 100
	double STranslationFactor2 =  5 ;
	double SRescaleFactor2 =  10 ;
	//int numberOfLevels[] = {1};
	//int numberOfBins[] = { 50 };
	for (int l = 0; l < sizeof(iterationsNumber2) / sizeof(iterationsNumber2[0]); l++)
	{
		for (int i = 0; i < sizeof(iterationsNumber1) / sizeof(iterationsNumber1[0]); i++)
		{
			for (int j = 0; j < sizeof(STranslationFactor) / sizeof(STranslationFactor[0]); j++)
			{
				for (int k = 0; k < sizeof(SRescaleFactor) / sizeof(SRescaleFactor[0]); k++)
				{					

					/*----------------------------------------------------
					            setting for 20151209 scans
					----------------------------------------------------*/
					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20151209/";
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20151209/";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20151209/MultiStage_base/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%sbase_ear_matched.nii", fixedfolder);
					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);

					/*----------------------------------------------------
					setting for 20160108  atlas-subject
					----------------------------------------------------*/
					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/MultiStage_second/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%ssecond_base_line_ear.nii", fixedfolder);
					////sprintf_s(fixedImageFile, "%sthird_line_injection_ear.nii", fixedfolder);
					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);

					/*----------------------------------------------------
					setting for 20160108 atlas-subject-subject
					----------------------------------------------------*/

					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/MultiStage_atlas13/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d_Similarity_IterationNum%d_TransScale%g_ScaleFactor%g", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l], iterationsNumber3, STranslationFactor2, SRescaleFactor2);
					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%sfirst_base_line_ear.nii.gz", fixedfolder);
					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);
					//char scan2ImageFile[300];
					//sprintf_s(scan2ImageFile, "%sthird_line_injection_ear.nii", fixedfolder);

					/*----------------------------------------------------
					setting for 20160307 scans
					----------------------------------------------------*/
					//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160307/";
					//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160307/";
					//char outputfolder[300];
					//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160307/MultiStage/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
					//char fixedImageFile[300];
					//sprintf_s(fixedImageFile, "%sRightEar.nii.gz", fixedfolder);
					//char movingImageFile[300];
					//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);

					/*----------------------------------------------------
					setting for 20160906 scans
					----------------------------------------------------*/
					
					char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";
					char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROI/";
					char outputfolder[300];
					sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/MultiStage/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%sminutes0.nii", fixedfolder);
					char movingImageFile[300];
					sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);
					
					//----------------------------------------------------------
					char parameter[50];
					sprintf_s(parameter, " %d %g %g %d %d %g %g", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l], iterationsNumber3, STranslationFactor2, SRescaleFactor2);

					char registeredImageFile[300];
					sprintf_s(registeredImageFile, " %sRegistrationOutput.nii", outputfolder);

					char differenceBeforeFile[300];
					sprintf_s(differenceBeforeFile, " %sdifferenceBeforeRegistration.nii", outputfolder);

					char differenceAfterFile[300];
					sprintf_s(differenceAfterFile, " %sdifferenceAfterRegistration.nii", outputfolder);

					char labelmapOutput_SL[300];
					sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_SM[300];
					sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_ST[300];
					sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", outputfolder);

					char labelmapOutput_SV[300];
					sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", outputfolder);

					char labelmapInput_SV[300];
					sprintf_s(labelmapInput_SV, " %satlasVolume_SV_labelmap.nii", inputfolder);
					char labelmapInput_ST[300];
					sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", inputfolder);
					char labelmapInput_SM[300];
					sprintf_s(labelmapInput_SM, " %satlasVolume_SM_labelmap.nii", inputfolder);
					char labelmapInput_SL[300];
					sprintf_s(labelmapInput_SL, " %satlasVolume_SL_labelmap.nii", inputfolder);

					//char labelmapInput_SV[300];
					//sprintf_s(labelmapInput_SV, " %sSV_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_ST[300];
					//sprintf_s(labelmapInput_ST, " %sST_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_SM[300];
					//sprintf_s(labelmapInput_SM, " %sSM_labelmap_registrationOutput.nii", inputfolder);
					//char labelmapInput_SL[300];
					//sprintf_s(labelmapInput_SL, " %sSL_labelmap_registrationOutput.nii", inputfolder);

				

					char command[3000] = "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/MultiStage2.exe ";

					strcat_s(command, fixedImageFile);
					strcat_s(command, movingImageFile);
					strcat_s(command, registeredImageFile);

					strcat_s(command, differenceBeforeFile);
					strcat_s(command, differenceAfterFile);

					strcat_s(command, parameter);
					//strcat_s(command, scan2ImageFile);
					strcat_s(command, labelmapInput_SV);
					strcat_s(command, labelmapOutput_SV);
					strcat_s(command, labelmapInput_ST);
					strcat_s(command, labelmapOutput_ST);
					strcat_s(command, labelmapInput_SM);
					strcat_s(command, labelmapOutput_SM);
					strcat_s(command, labelmapInput_SL);
					strcat_s(command, labelmapOutput_SL);

					std::cout << command << endl;
					std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

					std::system(command);
					std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

					std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;
				}
			}
		}
	}
	std::cout << "Done!" << std::endl;
	//getchar();
	return 0;
}