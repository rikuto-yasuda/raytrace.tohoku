////////////////////////////////////////////////////////////////////////
// tracer.inl
//	This is a library for software uses RAY-TRACING.
//	Copyright (C) 2005 Miyamoto Luisch
namespace rtc {

// ���C�g���[�V���O�����s���A����i�� //////////////////////////////////
double ray::take_a_step()
{
	// ��������x�N�g�������o���B
	vector& r  = m_rk.first;
	vector& k  = m_rk.second;

#ifndef NDEBUG
	intermediate i;
	update_intermediate(i,r,k);

	assert( i == m_im );
#endif

	const double dwdg = 1.0/calc_dGdw(m_im,r,k);
	const vector dgdr = calc_dGdr(m_im,r,k);
	const vector dgdk = calc_dGdk(m_im,r,k);

	m_drk.first  = -dgdk*dwdg; // -drdt
	m_drk.second =  dgdr*dwdg; //  dkdt

	// calc_dt() �̒��� update_intermediate()���Ăяo����
	// m_im ���X�V�����B
	const double dt = calc_dt(
		m_rk, m_drk,
		m_im
	);

	// ���ʂ��i�[���āA�����Ԃ��B
	m_drk.first  *= dt;
	m_drk.second *= dt;
	r += m_drk.first;
	k += m_drk.second;

	// �V�����_�ł̌��̏�Ԃ��`�F�b�N����B
	checkState(m_im,m_drk,r,k);

    if(reflectioncheck (r,m_drk)==1)
	{
		r -= m_drk.first;  //r,k = m_rk���ɖ߂�
		k -= m_drk.second;
		m_drk.first  /= dt;
		m_drk.second /= dt;

		const double dt = reflect_dt(m_rk,m_drk,m_im);

		m_drk.first  *= dt;
		m_drk.second *= dt;
		r += m_drk.first;
		k += m_drk.second;

		const vector n = reflect_n(m_rk,m_im);
		k = k-2*inner_prod(n,k)*n;
	}

	return dt; 

	
	//mainloop�ŌĂяo��������
	//dt=calc(), getDeltaR()=m_dkr.first, getDeltaK()=m_dkr.second
	//getR()=m_rk.first, getK()=m_rk.second �����ŎQ�Ƃ���Ă���̂Œl�̍X�V���\��
	//checkstate�Ōő̕��̔��f�����邽�߂����̒l��dt��Ԃ��O�ɕύX����K�v����
}

}// namespace rtc
