// wave_param.h: wave_parameter �N���X�̃C���^�[�t�F�C�X
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_WAVE_PARAM_H__C6B69AE0_831E_48CB_9621_49637B2FB295__INCLUDED_)
#define AFX_WAVE_PARAM_H__C6B69AE0_831E_48CB_9621_49637B2FB295__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

namespace rtc { // -----------------------------------------------------

	// wave_parameter --------------------------------------------------
	// �v�Z����g���̓������Ǘ�����N���X�B
	// �R���X�g���N�^�̈����ɁA�g���̃p�����[�^���w�肵�Đ�������B
	class wave_parameter  
	{
	public:
		// �萔
		enum wave_mode{
			LO_MODE = +1,
			RX_MODE = -1
		};

	public:
		wave_parameter(){};
		wave_parameter(
			enum wave_mode mode, // �g���̃��[�h�BLO_MODE��RX_MODE���w�肷��B
			double         freq, // �g���̎��g��[Hz]���w�肷��B
			double         precision,
			double         lstep,
			double         t_max,
			double         t_min
		);
		wave_parameter( const wave_parameter& r );

		virtual ~wave_parameter();

	public:
		//�l�̎Q��
		int       LO_or_RX()      const; // LO���[�h�Ȃ��+1�ARX���[�h�Ȃ�-1��Ԃ��B
		wave_mode getMode()       const; // enum wave_mode�l��Ԃ��B
		double    getStepLength() const; // ��X�e�b�v�Ői�ދ����́��Q�l�l����Ԃ��B
		double    getFreq()       const; // ���g����Ԃ��B
		double    getPrecision()  const; // ���x��Ԃ��B
		
		// ���Ԑ��x�͈̔͂�Ԃ��B
		const std::pair<double,double>& getTimeStep() const;

		// �l�̐ݒ�
		void setPrecision( const double prec );
		void setMode( wave_mode mode );
		void setStepLength( const double lstep );
		void setFreq( const double freq );

	private:
		wave_mode m_mode;
		double    m_freq;
		double    m_lstep;
		double    m_precision;
		std::pair<double,double> m_dt;
	};
}// namespace rtc; -----------------------------------------------------
#endif // !defined(AFX_WAVE_PARAM_H__C6B69AE0_831E_48CB_9621_49637B2FB295__INCLUDED_)
