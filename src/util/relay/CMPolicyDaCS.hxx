/*
	Definition of CMPolicyDaCS class

	Author: Benjamin Karl Bergen

	$Revision$
	$LastChangedBy$
	$LastChangedDate$

	vim: set ts=3 :
*/
#ifndef CMPolicyDaCS_hxx
#define CMPolicyDaCS_hxx

#include <iostream>
#include <assert.h>
#include <dacs.h>
#include <DaCSUtils.h>
#include <MPData.hxx>

/*!
	\struct CMPolicyDaCS CMPolicyDaCS.h
	\brief  provides...
*/
class CMPolicyDaCS
	{
	public:

		inline void init(int argc, char ** argv);
		inline void finalize();

		inline int global_id() { return id_; }
		inline int global_size() { return size_; }

		inline de_id_t peer_de() { return peer_de_; }
		inline dacs_process_id_t peer_pid() { return peer_pid_; }

	protected:

		CMPolicyDaCS()
			: id_(-1), size_(-1), peer_de_(DACS_DE_PARENT),
			peer_pid_(DACS_PID_PARENT) {}
		CMPolicyDaCS(const CMPolicyDaCS & mgr) {}
		~CMPolicyDaCS() {}

	private:

		int id_;
		int size_;

		DACS_ERR_T errcode_;
		de_id_t de_;
		de_id_t peer_de_;
		dacs_process_id_t peer_pid_;

	}; // class CMPolicyDaCS

void CMPolicyDaCS::init(int argc, char ** argv)
	{
		errcode_ = dacs_runtime_init(NULL, NULL, &de_);	
		process_dacs_errcode(errcode_);

		std::cerr << "accelerator process initialized" << std::endl;

		dacs_tag_t tag;
		errcode_ = dacs_tag_reserve(&tag, peer_de_, peer_pid_);
		process_dacs_errcode(errcode_);

		int data[2] __attribute__ ((aligned (16)));

		/*
		void * mem(NULL);
		posix_memalign((void **)(&mem), 16, 2*sizeof(int));
		int * data __attribute__ ((aligned (16))) = (int *)mem;
		*/

		errcode_ = dacs_recv(data, 2*sizeof(int), peer_de_,
			peer_pid_, tag, DACS_BYTE_SWAP_WORD);
		std::cerr << "accelerator called dacs_recv " << std::endl;
		process_dacs_errcode(errcode_);

		errcode_ = dacs_wait(tag);
		std::cerr << "accelerator called dacs_wait" << std::endl;
		process_dacs_errcode(errcode_);

		errcode_ = dacs_tag_release(&tag);
		std::cerr << "accelerator called dacs_tag_release" << std::endl;
		process_dacs_errcode(errcode_);

		id_ = data[0];
		size_ = data[1];

		std::cerr << "accel rank " << id_ << " size " << size_ << std::endl;

		//free(mem);
	} // CMPolicyDaCS::init

void CMPolicyDaCS::finalize()
	{
		errcode_ = dacs_runtime_exit();
		process_dacs_errcode(errcode_);
	} // CMPolicyDaCS::finalize

#endif // CMPolicyDaCS_hxx
