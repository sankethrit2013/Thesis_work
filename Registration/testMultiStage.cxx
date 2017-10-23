#include <iostream>
#include <string>
#include "stdio.h"
#include <chrono>

using namespace std;

int main(int argc, char *argv[])
{
	double STranslationFactor[] = {2};
	double SRescaleFactor[] = {5};
	double similarityScale[] = { 1 };
	//double SRescaleFactor[] = { 2, 5 };
	//int iterationsNumber1[] = {50, 200}; // 150
	int iterationsNumber1[] = { 150 };
	//int iterationsNumber2[] = {50, 100, 150}; // 0
	int iterationsNumber2[] = { 100 };
	//double similarityLearningRate[] = {0.01, 0.05, 0.1, 0.5};
	/*double affineLearningRate[] = { 0.005, 0.01, 0.05, 0.1, 0.5 };*/
	//double similarityLearningRate[] = { 0.01 };
	//double affineLearningRate[] = { 0.005};


	//int numberOfLevels[] = {1};
	//int numberOfBins[] = { 50 };
	for (int l = 0; l < sizeof(iterationsNumber2) / sizeof(iterationsNumber2[0]); l++)
	{
		for (int i = 0; i < sizeof(iterationsNumber1) / sizeof(iterationsNumber1[0]); i++)
		{
			for (int m = 0; m < sizeof(similarityScale) / sizeof(similarityScale[0]); m++)
			{
				/*for (int k = 0; k < sizeof(affineLearningRate) / sizeof(affineLearningRate[0]); k++)*/
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
						//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/MultiStage_first/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//char fixedImageFile[300];
						//sprintf_s(fixedImageFile, "%sfirst_base_line_ear.nii.gz", fixedfolder);
						////sprintf_s(fixedImageFile, "%ssecond_base_line_ear.nii", fixedfolder);
						////sprintf_s(fixedImageFile, "%sthird_line_injection_ear.nii", fixedfolder);
						//char movingImageFile[300];
						//sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);


						/*----------------------------------------------------
						setting for 20160108 subject-subject
						----------------------------------------------------*/
						//char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
						//char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/";
						//char outputfolder[300];
						//sprintf_s(outputfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160108/MultiStage_1to3/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//char fixedImageFile[300];
						//sprintf_s(fixedImageFile, "%sthird_line_injection_ear.nii", fixedfolder);
						//char movingImageFile[300];
						//sprintf_s(movingImageFile, "%sfirst_base_line_ear.nii.gz", inputfolder);


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
						setting for 20160719 scans
						----------------------------------------------------
						char* inputfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";
						char* fixedfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/";
						char* labelmapfolder = " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan1/labelmap_Affine/";
						char  labelmapfolder[400];
						sprintf_s(labelmapfolder, " C:/Users/zx4403/Dropbox/Research/ITK/ImageVolume/20160719/scan4/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						*/

						/*----------------------------------------------------
						setting for 20160906 scans
						----------------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/";
						//////std::cout << "Done inputfolder" << std::endl;

						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/WithROI/";
						//////std::cout << "Done fixedfolder" << std::endl;

						//////char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160719/";
						////////std::cout << "Done inputfolder" << std::endl;

						//////char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160719/";
						////////std::cout << "Done fixedfolder" << std::endl;

						//////char outputfolder[300];
						//////sprintf_s(outputfolder, " N:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/MultiStage_other/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//////char fixedImageFile[300];
						//////sprintf_s(fixedImageFile, "%sbase_line.nii", fixedfolder);
						//////char movingImageFile[300];
						//////sprintf_s(movingImageFile, "%satlasVolume_unlabeled.nii", inputfolder);
						//////char* labelmapfolder = " N:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/baseline/labelmap_Affine/";
						//char  labelmapfolder[400];
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						/////*char  labelmapfolder[400];
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/old_testing/labelmap_Affine/baseline_Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);*/

						//////----------------------------------------------------------
						////
						////
						/////*char outputfolder[300];
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/old_testing/MultiStage/baseline_Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);*/

						//char outputfolder[300];
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160906/MultiStage/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						/*-------------------------------------------------------------------------------------------------*/

						/*--------------------------------------------
						20160629 scans
						---------------------------------------------*/
						/*char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/";
						char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/scan1/";
						char  labelmapfolder[400];
						sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/labelmap_Affine/atlas_Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						char outputfolder[300];
						sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20160629/MultiStage/atlas_Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);*/

						/*-------------------------------------------*/

						///*-------------------------------------------
						//20170111 scans
						//---------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/";
						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/scan1/";

						//char  labelmapfolder[400];
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						//char outputfolder[300];
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170111/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFactor%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						///*-------------------------------------------*/

						///*-------------------------------------------
						//20170123 scans
						//---------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/";
						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/scan1/";

						//char  labelmapfolder[400];
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/labelmap_Affine/Similarity_IterationNum%d_SimLR%g_AffineLR%g_Affine_IterationNum%d", iterationsNumber1[i], similarityLearningRate[j], affineLearningRate[k], iterationsNumber2[l]);
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/labelmap_Affine/CL01638_Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						//char outputfolder[300];
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/Multistage/Similarity_IterationNum%d_SimLR%g_AffineLR%g_Affine_IterationNum%d", iterationsNumber1[i], similarityLearningRate[j], affineLearningRate[k], iterationsNumber2[l]);
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170123/Multistage/CL01638_Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						///*-------------------------------------------*/

						///*-------------------------------------------
						//20170606 scans
						//---------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";
						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/";

						//char  labelmapfolder[400];
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/labelmap_Affine/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						//char outputfolder[300];
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Multistage/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						///*-------------------------------------------*/

						///*-------------------------------------------
						//20170612 scans
						//---------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/";
						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/scan1/";

						//char  labelmapfolder[400];
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/labelmap_Affine/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						//char outputfolder[300];
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170612/Multistage/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						///*-------------------------------------------*/

						///*-------------------------------------------
						//20170619 scans
						//---------------------------------------------*/
						//char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/";
						//char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/scan1/";

						//char  labelmapfolder[400];
						////sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/labelmap_Affine/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						//char outputfolder[300];
						////sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Multistage/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						///*-------------------------------------------*/

						/*-------------------------------------------*/

						/*-------------------------------------------
						20171017 scans
						---------------------------------------------*/
						char* inputfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/";
						char* fixedfolder = " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/scan1/";

						char  labelmapfolder[400];
						//sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/labelmap_Affine/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						sprintf_s(labelmapfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/labelmap_Affine/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						char outputfolder[300];
						//sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20170606/Multistage/Similarity_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);
						sprintf_s(outputfolder, " E:/Sanketh_Moudgalya_Directory/Nate/ImageVolume/20171017/Multistage/Similarity_Scale%g_IterationNum%d_TransScale%g_ScaleFact%g_Affine_IterationNum%d", similarityScale[m], iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l]);

						/*-------------------------------------------*/


						char fixedImageFile[300];
						//sprintf_s(fixedImageFile, "%sscan1rescale.nii", inputfolder);
						sprintf_s(fixedImageFile, "%sscan1.nii", fixedfolder);

						char movingImageFile[300];
						//sprintf_s(movingImageFile, "%sscan1atlas.nii", fixedfolder);
						sprintf_s(movingImageFile, "%satlas_bspline3.nii", inputfolder);
						//sprintf_s(movingImageFile, "%sscan1.nii", fixedfolder);

						char parameter[50];
						sprintf_s(parameter, " %d %g %g %d %g ", iterationsNumber1[i], STranslationFactor[j], SRescaleFactor[k], iterationsNumber2[l], similarityScale[m]);
						//sprintf_s(parameter, " %d %g %g %g %g %d", iterationsNumber1[i], STranslationFactor, SRescaleFactor, similarityLearningRate[j], affineLearningRate[k], iterationsNumber2[l]);

						char registeredImageFile[300];
						sprintf_s(registeredImageFile, " %sRegistrationOutput.nii", outputfolder);

						char differenceBeforeFile[300];
						sprintf_s(differenceBeforeFile, " %sdifferenceBeforeRegistration.nii", outputfolder);

						char differenceAfterFile[300];
						sprintf_s(differenceAfterFile, " %sdifferenceAfterRegistration.nii.gz", outputfolder);

#pragma region labelMapFileName
						//label map output
						char labelmapOutput_SL[300];
						sprintf_s(labelmapOutput_SL, " %sSL_labelmap_registrationOutput.nii", labelmapfolder);

						char labelmapOutput_SM[300];
						sprintf_s(labelmapOutput_SM, " %sSM_labelmap_registrationOutput.nii", labelmapfolder);
						/*char labelmapOutput_SM_e1[300];
						sprintf_s(labelmapOutput_SM_e1, " %sSM_e1_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SM_e2[300];
						sprintf_s(labelmapOutput_SM_e2, " %sSM_e2_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SM_e3[300];
						sprintf_s(labelmapOutput_SM_e3, " %sSM_e3_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SM_e4[300];
						sprintf_s(labelmapOutput_SM_e4, " %sSM_e4_labelmap_registrationOutput.nii", labelmapfolder);*/

						char labelmapOutput_ST[300];
						sprintf_s(labelmapOutput_ST, " %sST_labelmap_registrationOutput.nii", labelmapfolder);
						/*char labelmapOutput_ST_e1[300];
						sprintf_s(labelmapOutput_ST_e1, " %sST_e1_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_ST_e2[300];
						sprintf_s(labelmapOutput_ST_e2, " %sST_e2_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_ST_e3[300];
						sprintf_s(labelmapOutput_ST_e3, " %sST_e3_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_ST_e4[300];
						sprintf_s(labelmapOutput_ST_e4, " %sST_e4_labelmap_registrationOutput.nii", labelmapfolder);*/

						char labelmapOutput_SV[300];
						sprintf_s(labelmapOutput_SV, " %sSV_labelmap_registrationOutput.nii", labelmapfolder);
						/*char labelmapOutput_SV_e1[300];
						sprintf_s(labelmapOutput_SV_e1, " %sSV_e1_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SV_e2[300];
						sprintf_s(labelmapOutput_SV_e2, " %sSV_e2_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SV_e3[300];
						sprintf_s(labelmapOutput_SV_e3, " %sSV_e3_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SV_e4[300];
						sprintf_s(labelmapOutput_SV_e4, " %sSV_e4_labelmap_registrationOutput.nii", labelmapfolder);
						char labelmapOutput_SV_e5[300];
						sprintf_s(labelmapOutput_SV_e5, " %sSV_e5_labelmap_registrationOutput.nii", labelmapfolder);*/
						char labelmapOutput_CA[300];
						sprintf_s(labelmapOutput_CA, " %sCA_labelmap_registrationOutput.nii", labelmapfolder);

						/*char labelmapOutput_Comb[300];
						sprintf_s(labelmapOutput_Comb, " %slabelmap_CL_combined01638_registrationOutput.nii", labelmapfolder);*/

						//label map input
						char labelmapInput_SV[300];
						sprintf_s(labelmapInput_SV, " %satlasVolume_SV_labelmap.nii", inputfolder);
						/*char labelmapInput_SV_e1[300];
						sprintf_s(labelmapInput_SV_e1, " %satlasVolume_SV_labelmap_e1.nii", inputfolder);
						char labelmapInput_SV_e2[300];
						sprintf_s(labelmapInput_SV_e2, " %satlasVolume_SV_labelmap_e2.nii", inputfolder);
						char labelmapInput_SV_e3[300];
						sprintf_s(labelmapInput_SV_e3, " %satlasVolume_SV_labelmap_e3.nii", inputfolder);
						char labelmapInput_SV_e4[300];
						sprintf_s(labelmapInput_SV_e4, " %satlasVolume_SV_labelmap_e4.nii", inputfolder);
						char labelmapInput_SV_e5[300];
						sprintf_s(labelmapInput_SV_e5, " %satlasVolume_SV_labelmap_e5.nii", inputfolder);*/

						char labelmapInput_ST[300];
						sprintf_s(labelmapInput_ST, " %satlasVolume_ST_labelmap.nii", inputfolder);
						/*char labelmapInput_ST_e1[300];
						sprintf_s(labelmapInput_ST_e1, " %satlasVolume_ST_labelmap_e1.nii", inputfolder);
						char labelmapInput_ST_e2[300];
						sprintf_s(labelmapInput_ST_e2, " %satlasVolume_ST_labelmap_e2.nii", inputfolder);
						char labelmapInput_ST_e3[300];
						sprintf_s(labelmapInput_ST_e3, " %satlasVolume_ST_labelmap_e3.nii", inputfolder);
						char labelmapInput_ST_e4[300];
						sprintf_s(labelmapInput_ST_e4, " %satlasVolume_ST_labelmap_e4.nii", inputfolder);*/

						char labelmapInput_SM[300];
						sprintf_s(labelmapInput_SM, " %satlasVolume_SM_labelmap.nii", inputfolder);
						/*char labelmapInput_SM_e1[300];
						sprintf_s(labelmapInput_SM_e1, " %satlasVolume_SM_labelmap_e1.nii", inputfolder);
						char labelmapInput_SM_e2[300];
						sprintf_s(labelmapInput_SM_e2, " %satlasVolume_SM_labelmap_e2.nii", inputfolder);
						char labelmapInput_SM_e3[300];
						sprintf_s(labelmapInput_SM_e3, " %satlasVolume_SM_labelmap_e3.nii", inputfolder);
						char labelmapInput_SM_e4[300];
						sprintf_s(labelmapInput_SM_e4, " %satlasVolume_SM_labelmap_e4.nii", inputfolder);*/

						char labelmapInput_SL[300];
						sprintf_s(labelmapInput_SL, " %satlasVolume_SL_labelmap.nii", inputfolder);
						char labelmapInput_CA[300];
						sprintf_s(labelmapInput_CA, " %sCA_labelmap.nii", inputfolder);

						/*char labelmapInput_Comb[300];
						sprintf_s(labelmapInput_Comb, " %slabelmap_combined01638.nii", inputfolder);*/
#pragma endregion



						//char labelmapInput_SV[300];
						//sprintf_s(labelmapInput_SV, " %sSV_labelmap_registrationOutput.nii", inputfolder);
						//char labelmapInput_ST[300];
						//sprintf_s(labelmapInput_ST, " %sST_labelmap_registrationOutput.nii", inputfolder);
						//char labelmapInput_SM[300];
						//sprintf_s(labelmapInput_SM, " %sSM_labelmap_registrationOutput.nii", inputfolder);
						//char labelmapInput_SL[300];
						//sprintf_s(labelmapInput_SL, " %sSL_labelmap_registrationOutput.nii", inputfolder);



						char command[8000] = "";

						strcat_s(command, "E:/Sanketh_Moudgalya_Directory/Nate/C++_codes/built_codes/Regi/Debug/MultiStage.exe ");

						strcat_s(command, fixedImageFile);
						strcat_s(command, movingImageFile);
						strcat_s(command, registeredImageFile);

						strcat_s(command, differenceBeforeFile);

						//std::cout << "user command" << std::endl;
						//std::cout << command << std::endl;
						strcat_s(command, differenceAfterFile);
						strcat_s(command, parameter);
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

						/*strcat_s(command, labelmapInput_SV_e1);
						strcat_s(command, labelmapOutput_SV_e1);
						strcat_s(command, labelmapInput_SV_e2);
						strcat_s(command, labelmapOutput_SV_e2);
						strcat_s(command, labelmapInput_SV_e3);
						strcat_s(command, labelmapOutput_SV_e3);
						strcat_s(command, labelmapInput_SV_e4);
						strcat_s(command, labelmapOutput_SV_e4);
						strcat_s(command, labelmapInput_SV_e5);
						strcat_s(command, labelmapOutput_SV_e5);

						strcat_s(command, labelmapInput_ST_e1);
						strcat_s(command, labelmapOutput_ST_e1);
						strcat_s(command, labelmapInput_ST_e2);
						strcat_s(command, labelmapOutput_ST_e2);
						strcat_s(command, labelmapInput_ST_e3);
						strcat_s(command, labelmapOutput_ST_e3);
						strcat_s(command, labelmapInput_ST_e4);
						strcat_s(command, labelmapOutput_ST_e4);

						strcat_s(command, labelmapInput_SM_e1);
						strcat_s(command, labelmapOutput_SM_e1);
						strcat_s(command, labelmapInput_SM_e2);
						strcat_s(command, labelmapOutput_SM_e2);
						strcat_s(command, labelmapInput_SM_e3);
						strcat_s(command, labelmapOutput_SM_e3);
						strcat_s(command, labelmapInput_SM_e4);
						strcat_s(command, labelmapOutput_SM_e4);*/

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
	}
	return 0;
}