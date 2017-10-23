#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>
#include <itkMacro.h>

using namespace std;

int main(int argc, char *argv[])
{
	double STranslationFactor[] = { 2 };
	double SRescaleFactor[] = { 5 };
	int iterationsNumber1[] = { 50 };
	int iterationsNumber2[] = { 20 };

	for (int l = 0; l < sizeof(iterationsNumber2) / sizeof(iterationsNumber2[0]); l++)
	{
		for (int i = 0; i < sizeof(iterationsNumber1) / sizeof(iterationsNumber1[0]); i++)
		{
			/*for (int j = 0; j < sizeof(similarityLearningRate) / sizeof(similarityLearningRate[0]); j++)
			{
			for (int k = 0; k < sizeof(affineLearningRate) / sizeof(affineLearningRate[0]); k++)*/
			for (int j = 0; j < sizeof(STranslationFactor) / sizeof(STranslationFactor[0]); j++)
			{
				for (int k = 0; k < sizeof(SRescaleFactor) / sizeof(SRescaleFactor[0]); k++)
				{
					/*-------------------------------------------
					20170123 scans
					---------------------------------------------*/
					char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/";
					char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/";

					char  labelmapfolder[400];
					sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/labelmap_Affine/Resampled_Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

					char outputfolder[300];
					sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/Multistage/New2D_Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

					char fixedImageFile[300];
					sprintf_s(fixedImageFile, "%suct1.tif", inputfolder);

					char movingImageFile[300];
					sprintf_s(movingImageFile, "%satlas1.tif", inputfolder);

					char parameter[50];
					sprintf_s(parameter, " %d %g %g %d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

					char registeredImageFile[300];
					sprintf_s(registeredImageFile, " %sRegistrationOutput.jpg", outputfolder);

					char differenceBeforeFile[300];
					sprintf_s(differenceBeforeFile, " %sdifferenceBeforeRegistration.nii", outputfolder);

					char differenceAfterFile[300];
					sprintf_s(differenceAfterFile, " %sdifferenceAfterRegistration.jpg", outputfolder);

#pragma region labelMapFileName

					// labelmap output
					char labelmapOutput_ST[300];
					sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", labelmapfolder);
					// labelmap input
					char labelmapInput_ST[300];
					sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", inputfolder);

#pragma endregion

					char command[8000] = "";

					strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/MasoumehCode/Debug/Multistage2D.exe ");

					strcat_s(command, fixedImageFile);
					strcat_s(command, movingImageFile);
					strcat_s(command, registeredImageFile);

					strcat_s(command, differenceBeforeFile);
					strcat_s(command, differenceAfterFile);
					strcat_s(command, parameter);
					strcat_s(command, labelmapInput_ST);
					strcat_s(command, labelmapOutput_ST);

					std::cout << "Creating command" << std::endl;
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
//int main()
//{
//	double translationFactor[] = { 0.1 };
//	int iterationsNumber[] = { 100 };
//	double learningRate[] = { 0.1 };
//	double minimumStepLength[] = { 0.00001 };
//	for (size_t t = 1; t <= 6; t++)
//	{
//		if (t != 1)
//		{
//			//sanketh code
//			std::cout << "scan" << t << " computing" << std::endl;
//
//			for (int l = 0; l < sizeof(minimumStepLength) / sizeof(minimumStepLength[0]); l++)
//			{
//				for (int j = 0; j < sizeof(learningRate) / sizeof(learningRate[0]); j++)
//				{
//					for (int k = 0; k < sizeof(iterationsNumber) / sizeof(iterationsNumber[0]); k++)
//					{
//						for (int i = 0; i < sizeof(translationFactor) / sizeof(translationFactor[0]); i++)
//						{
//							char* scansfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/";
//							char* baselinefolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/";
//							char otherscansfolder[300];
//							sprintf_s(otherscansfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan%d/", t);
//
//
//							char* labelmapInputFolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/labelmap_Deformable_Grid10_Iter150/";
//							char labelmapOutputFolder[300];
//							sprintf_s(labelmapOutputFolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan%d/labelmap_Deformable_Grid10_Iter150_fromscan1/", t);
//							/*-----------------------------------------------------------------------------------------------------------------*/
//
//							char fixedImageFile[300];
//							sprintf_s(fixedImageFile, "%sscan%d.nii", otherscansfolder, t);
//							char movingImageFile[300];
//							sprintf_s(movingImageFile, "%sscan1.nii", baselinefolder, t);
//
//							char ResampleInputFile[300];
//							sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension10_Iterations150RegistrationOutput.nii", scansfolder);
//							//sprintf_s(ResampleInputFile, "%sscan1/BSpline_GridNodeinOneDimension10_Iterations100RegistrationOutput.nii", scansfolder);
//							char ResampleOutputFile[300];
//							sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations150RegistrationOutput.nii", scansfolder, t);
//							//sprintf_s(ResampleOutputFile, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations100_RegistrationOutput.nii", scansfolder, t);
//							//sprintf_s(ResampleOutputFile, "%sRegistrationOutput.nii", scansfolder);
//
//
//							char differenceAfterFile_scan2scan[300];
//							sprintf_s(differenceAfterFile_scan2scan, "%sscan%d/Rigid_toScan1_iterationNum%d_translationFactor%d_learningRate%g_minStepLength%g_differenceAfterRegistration.nii", scansfolder, t, iterationsNumber[k], translationFactor[i], learningRate[j], minimumStepLength[l]);
//
//							char differenceAfterFile_scan2atlas[300];
//							/*sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension8_Iterations50_differenceAfterRegistration.nii", scansfolder, t);*/
//							sprintf_s(differenceAfterFile_scan2atlas, "%sscan%d/BSpline_fromScan1_GridNodeinOneDimension10_Iterations150_differenceAfterRegistration.nii", scansfolder, t);
//
//							char parameter[50];
//							sprintf_s(parameter, " %g %d %g %g", translationFactor[i], iterationsNumber[k], learningRate[j], minimumStepLength[l]);
//
//							char command[8000] = "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/RegistrationRigid3D.exe ";
//
//							strcat_s(command, fixedImageFile);
//							strcat_s(command, movingImageFile);
//							strcat_s(command, ResampleInputFile);
//
//							strcat_s(command, ResampleOutputFile);
//							strcat_s(command, differenceAfterFile_scan2scan);
//							strcat_s(command, differenceAfterFile_scan2atlas);
//							strcat_s(command, parameter);
//
//							cout << command << endl;
//							std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//
//							system(command);
//							std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//
//							std::cout << "Running time = " << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count() << " minutes" << std::endl;
//
//						}
//					}
//				}
//			}
//		}
//	}
//	return 0;
//}