
        const __m512i gather_indices = _mm512_load_epi32( gather_indices_arr );

        // Build mask:
        const __mmask16 gather_mask = _mm512_cmplt_epi32_mask(linear, _mm512_set1_epi32(to_process));

        // Perform gathers:
        { // ex {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ex, 4);
            _mm512_mask_store_ps(&ex[0], gather_mask, gathered_values);
        } // }}}
        { // dexdy {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dexdy, 4);
            _mm512_mask_store_ps(&dexdy[0], gather_mask, gathered_values);
        } // }}}
        { // dexdz {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dexdz, 4);
            _mm512_mask_store_ps(&dexdz[0], gather_mask, gathered_values);
        } // }}}
        { // d2exdydz {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2exdydz, 4);
            _mm512_mask_store_ps(&d2exdydz[0], gather_mask, gathered_values);
        } // }}}

        { // ey {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ey, 4);
            _mm512_mask_store_ps(&ey[0], gather_mask, gathered_values);
        } // }}}
        { // deydz {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].deydz, 4);
            _mm512_mask_store_ps(&deydz[0], gather_mask, gathered_values);
        } // }}}
        { // deydx {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].deydx, 4);
            _mm512_mask_store_ps(&deydx[0], gather_mask, gathered_values);
        } // }}}
        { // d2eydzdx {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2eydzdx, 4);
            _mm512_mask_store_ps(&d2eydzdx[0], gather_mask, gathered_values);
        } // }}}

        { // ez {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].ez, 4);
            _mm512_mask_store_ps(&ez[0], gather_mask, gathered_values);
        } // }}}
        { // dezdx {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dezdx, 4);
            _mm512_mask_store_ps(&dezdx[0], gather_mask, gathered_values);
        } // }}}
        { // dezdy {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dezdy, 4);
            _mm512_mask_store_ps(&dezdy[0], gather_mask, gathered_values);
        } // }}}
        { // d2ezdxdy {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].d2ezdxdy, 4);
            _mm512_mask_store_ps(&d2ezdxdy[0], gather_mask, gathered_values);
        } // }}}

        { // cbx {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cbx, 4);
            _mm512_mask_store_ps(&_cbx[0], gather_mask, gathered_values);
        } // }}}
        { // dcbxdx {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbxdx, 4);
            _mm512_mask_store_ps(&dcbxdx[0], gather_mask, gathered_values);
        } // }}}
        { // cby {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cby, 4);
            _mm512_mask_store_ps(&_cby[0], gather_mask, gathered_values);
        } // }}}
        { // dcbydy {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbydy, 4);
            _mm512_mask_store_ps(&dcbydy[0], gather_mask, gathered_values);
        } // }}}

        { // cbz {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].cbz, 4);
            _mm512_mask_store_ps(&_cbz[0], gather_mask, gathered_values);
        } // }}}
        { // dcbzdz {{{
            const __m512 gathered_values = _mm512_mask_i32gather_ps(zeroes, gather_mask, gather_indices, &f0[0].dcbzdz, 4);
            _mm512_mask_store_ps(&dcbzdz[0], gather_mask, gathered_values);
        } // }}}
        // }}}
        // */
        }
