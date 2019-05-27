// Copyright: (2012-2015) Ben Strasser <code@ben-strasser.net>
// License: BSD-3
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef CSV_H
#define CSV_H

#include <cstring>
#include <utility>
#include <cstdio>
#include <exception>
#include <memory>
#include <cerrno>
#include <istream>

/*
 * rapidcsv.h
 *
 * URL:      https://github.com/d99kris/rapidcsv
 * Version:  4.0
 *
 * Copyright (C) 2017-2019 Kristofer Berggren
 * All rights reserved.
 *
 * rapidcsv is distributed under the BSD 3-Clause license, see LICENSE for details.
 *
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#ifdef HAS_CODECVT
#include <codecvt>
#endif
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#if defined(_MSC_VER)
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

namespace rapidcsv
{
#if defined(_MSC_VER)
	static const bool sPlatformHasCR = true;
#else
	static const bool sPlatformHasCR = false;
#endif

	/**
	 * @brief     Datastructure holding parameters controlling how invalid numbers (including
	 *            empty strings) should be handled.
	 */
	struct ConverterParams
	{
		/**
		 * @brief   Constructor
		 * @param   pHasDefaultConverter  specifies if conversion of non-numerical strings shall be
		 *                                converted to a default numerical value, instead of causing
		 *                                an exception to be thrown (default).
		 * @param   pDefaultFloat         floating-point default value to represent invalid numbers.
		 * @param   pDefaultInteger       integer default value to represent invalid numbers.
		 */
		explicit ConverterParams(const bool pHasDefaultConverter = false,
			const long double pDefaultFloat = std::numeric_limits<long double>::signaling_NaN(),
			const long long pDefaultInteger = 0)
			: mHasDefaultConverter(pHasDefaultConverter)
			, mDefaultFloat(pDefaultFloat)
			, mDefaultInteger(pDefaultInteger)
		{
		}

		/**
		 * @brief   specifies if conversion of non-numerical strings shall be converted to a default
		 *          numerical value, instead of causing an exception to be thrown (default).
		 */
		bool mHasDefaultConverter;

		/**
		 * @brief   floating-point default value to represent invalid numbers.
		 */
		long double mDefaultFloat;

		/**
		 * @brief   integer default value to represent invalid numbers.
		 */
		long long mDefaultInteger;
	};

	/**
	 * @brief     Exception thrown when attempting to access Document data in a datatype which
	 *            is not supported by the Converter class.
	 */
	class no_converter : public std::exception
	{
		/**
		 * @brief   Provides details about the exception
		 * @returns an explanatory string
		 */
		virtual const char* what() const throw()
		{
			return "unsupported conversion datatype";
		}
	};

	/**
	 * @brief     Class providing conversion to/from numerical datatypes and strings. Only
	 *            intended for rapidcsv internal usage, but exposed externally to allow
	 *            specialization for custom datatype conversions.
	 */
	template<typename T>
	class Converter
	{
	public:
		/**
		 * @brief   Constructor
		 * @param   pConverterParams      specifies how conversion of non-numerical values to
		 *                                numerical datatype shall be handled.
		 */
		Converter(const ConverterParams& pConverterParams)
			: mConverterParams(pConverterParams)
		{
		}

		/**
		 * @brief   Converts numerical value to string representation.
		 * @param   pVal                  numerical value
		 * @param   pStr                  output string
		 */
		void ToStr(const T& pVal, std::string& pStr) const
		{
			if (typeid(T) == typeid(int) ||
				typeid(T) == typeid(long) ||
				typeid(T) == typeid(long long) ||
				typeid(T) == typeid(unsigned) ||
				typeid(T) == typeid(unsigned long) ||
				typeid(T) == typeid(unsigned long long) ||
				typeid(T) == typeid(float) ||
				typeid(T) == typeid(double) ||
				typeid(T) == typeid(long double) ||
				typeid(T) == typeid(char))
			{
				std::ostringstream out;
				out << pVal;
				pStr = out.str();
			}
			else
			{
				throw no_converter();
			}
		}

		/**
		 * @brief   Converts string holding a numerical value to numerical datatype representation.
		 * @param   pVal                  numerical value
		 * @param   pStr                  output string
		 */
		void ToVal(const std::string& pStr, T& pVal) const
		{
			try
			{
				if (typeid(T) == typeid(int))
				{
					pVal = static_cast<T>(std::stoi(pStr));
					return;
				}
				else if (typeid(T) == typeid(long))
				{
					pVal = static_cast<T>(std::stol(pStr));
					return;
				}
				else if (typeid(T) == typeid(long long))
				{
					pVal = static_cast<T>(std::stoll(pStr));
					return;
				}
				else if (typeid(T) == typeid(unsigned))
				{
					pVal = static_cast<T>(std::stoul(pStr));
					return;
				}
				else if (typeid(T) == typeid(unsigned long))
				{
					pVal = static_cast<T>(std::stoul(pStr));
					return;
				}
				else if (typeid(T) == typeid(unsigned long long))
				{
					pVal = static_cast<T>(std::stoull(pStr));
					return;
				}
			}
			catch (...)
			{
				if (!mConverterParams.mHasDefaultConverter)
				{
					throw;
				}
				else
				{
					pVal = static_cast<T>(mConverterParams.mDefaultInteger);
					return;
				}
			}

			try
			{
				if (typeid(T) == typeid(float))
				{
					pVal = static_cast<T>(std::stof(pStr));
					return;
				}
				else if (typeid(T) == typeid(double))
				{
					pVal = static_cast<T>(std::stod(pStr));
					return;
				}
				else if (typeid(T) == typeid(long double))
				{
					pVal = static_cast<T>(std::stold(pStr));
					return;
				}
			}
			catch (...)
			{
				if (!mConverterParams.mHasDefaultConverter)
				{
					throw;
				}
				else
				{
					pVal = static_cast<T>(mConverterParams.mDefaultFloat);
					return;
				}
			}

			if (typeid(T) == typeid(char))
			{
				pVal = static_cast<T>(pStr[0]);
				return;
			}
			else
			{
				throw no_converter();
			}
		}

	private:
		const ConverterParams& mConverterParams;
	};

	/**
	 * @brief     Specialized implementation handling string to string conversion.
	 * @param     pVal                  string
	 * @param     pStr                  string
	 */
	template<>
	inline void Converter<std::string>::ToStr(const std::string& pVal, std::string& pStr) const
	{
		pStr = pVal;
	}

	/**
	 * @brief     Specialized implementation handling string to string conversion.
	 * @param     pVal                  string
	 * @param     pStr                  string
	 */
	template<>
	inline void Converter<std::string>::ToVal(const std::string& pStr, std::string& pVal) const
	{
		pVal = pStr;
	}

	/**
	 * @brief     Datastructure holding parameters controlling which row and column should be
	 *            treated as labels.
	 */
	struct LabelParams
	{
		/**
		 * @brief   Constructor
		 * @param   pColumnNameIdx        specifies the zero-based row index of the column labels, setting
		 *                                it to -1 prevents column lookup by label name, and gives access
		 *                                to all rows as document data.
		 * @param   pRowNameIdx           specifies the zero-based column index of the row labels, setting
		 *                                it to -1 prevents row lookup by label name, and gives access
		 *                                to all columns as document data.
		 */
		explicit LabelParams(const int pColumnNameIdx = 0, const int pRowNameIdx = 0)
			: mColumnNameIdx(pColumnNameIdx)
			, mRowNameIdx(pRowNameIdx)
		{
		}

		/**
		 * @brief   specifies the zero-based row index of the column labels.
		 */
		int mColumnNameIdx;

		/**
		 * @brief   specifies the zero-based column index of the row labels.
		 */
		int mRowNameIdx;
	};

	/**
	 * @brief     Datastructure holding parameters controlling how the CSV data fields are separated.
	 */
	struct SeparatorParams
	{
		/**
		 * @brief   Constructor
		 * @param   pSeparator            specifies the column separator (default ',').
		 * @param   pHasCR                specifies whether a new document (i.e. not an existing document read)
		 *                                should use CR/LF instead of only LF (default is to use standard
		 *                                behavior of underlying platforms - CR/LF for Win, and LF for others).
		 */
		explicit SeparatorParams(const char pSeparator = ',', const bool pHasCR = sPlatformHasCR)
			: mSeparator(pSeparator)
			, mHasCR(pHasCR)
		{
		}

		/**
		 * @brief   specifies the column separator.
		 */
		char mSeparator;

		/**
		 * @brief   specifies whether new documents should use CR/LF instead of LF.
		 */
		bool mHasCR;
	};

	/**
	 * @brief     Class representing a CSV document.
	 */
	class Document
	{
	public:
		/**
		 * @brief   Constructor
		 * @param   pPath                 specifies the path of an existing CSV-file to populate the Document
		 *                                data with.
		 * @param   pLabelParams          specifies which row and column should be treated as labels.
		 * @param   pSeparatorParams      specifies which field and row separators should be used.
		 * @param   pConverterParams      specifies how invalid numbers (including empty strings) should be
		 *                                handled.
		 */
		explicit Document(const std::string& pPath = std::string(),
			const LabelParams& pLabelParams = LabelParams(),
			const SeparatorParams& pSeparatorParams = SeparatorParams(),
			const ConverterParams& pConverterParams = ConverterParams())
			: mPath(pPath)
			, mLabelParams(pLabelParams)
			, mSeparatorParams(pSeparatorParams)
			, mConverterParams(pConverterParams)
		{
			if (!mPath.empty())
			{
				ReadCsv();
			}
		}

		/**
		 * @brief   Constructor
		 * @param   pStream               specifies an input stream to read CSV data from.
		 * @param   pLabelParams          specifies which row and column should be treated as labels.
		 * @param   pSeparatorParams      specifies which field and row separators should be used.
		 * @param   pConverterParams      specifies how invalid numbers (including empty strings) should be
		 *                                handled.
		 */
		explicit Document(std::istream& pStream,
			const LabelParams& pLabelParams = LabelParams(),
			const SeparatorParams& pSeparatorParams = SeparatorParams(),
			const ConverterParams& pConverterParams = ConverterParams())
			: mPath()
			, mLabelParams(pLabelParams)
			, mSeparatorParams(pSeparatorParams)
			, mConverterParams(pConverterParams)
		{
			ReadCsv(pStream);
		}


		/**
		 * @brief   Copy constructor
		 * @param   pDocument             specifies the Document instance to copy.
		 */
		explicit Document(const Document& pDocument)
			: mPath(pDocument.mPath)
			, mLabelParams(pDocument.mLabelParams)
			, mSeparatorParams(pDocument.mSeparatorParams)
			, mConverterParams(pDocument.mConverterParams)
			, mData(pDocument.mData)
			, mColumnNames(pDocument.mColumnNames)
			, mRowNames(pDocument.mRowNames)
		{
		}

		/**
		 * @brief   Read Document data from file.
		 * @param   pPath                 specifies the path of an existing CSV-file to populate the Document
		 *                                data with.
		 */
		void Load(const std::string& pPath)
		{
			mPath = pPath;
			ReadCsv();
		}

		/**
		 * @brief   Write Document data to file.
		 * @param   pPath                 optionally specifies the path where the CSV-file will be created
		 *                                (if not specified, the original path provided when creating or
		 *                                loading the Document data will be used).
		 */
		void Save(const std::string& pPath = std::string())
		{
			if (!pPath.empty())
			{
				mPath = pPath;
			}
			WriteCsv();
		}

		/**
		 * @brief   Write Document data to stream.
		 * @param   pStream               specifies an output stream to write the data to.
		 */
		void Save(std::ostream& pStream)
		{
			WriteCsv(pStream);
		}

		/**
		 * @brief   Get column by index.
		 * @param   pColumnIdx            zero-based column index.
		 * @returns vector of column data.
		 */
		template<typename T>
		std::vector<T> GetColumn(const size_t pColumnIdx) const
		{
			const ssize_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			std::vector<T> column;
			Converter<T> converter(mConverterParams);
			for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
			{
				if (std::distance(mData.begin(), itRow) > mLabelParams.mColumnNameIdx)
				{
					T val;
					converter.ToVal(itRow->at(columnIdx), val);
					column.push_back(val);
				}
			}
			return column;
		}

		/**
		 * @brief   Get column by name.
		 * @param   pColumnName           column label name.
		 * @returns vector of column data.
		 */
		template<typename T>
		std::vector<T> GetColumn(const std::string& pColumnName) const
		{
			const ssize_t columnIdx = GetColumnIdx(pColumnName);
			if (columnIdx < 0)
			{
				throw std::out_of_range("column not found: " + pColumnName);
			}
			return GetColumn<T>(columnIdx);
		}

		/**
		 * @brief   Set column by index.
		 * @param   pColumnIdx            zero-based column index.
		 * @param   pColumn               vector of column data.
		 */
		template<typename T>
		void SetColumn(const size_t pColumnIdx, const std::vector<T>& pColumn)
		{
			const size_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);

			while (pColumn.size() + (mLabelParams.mColumnNameIdx + 1) > GetDataRowCount())
			{
				std::vector<std::string> row;
				row.resize(GetDataColumnCount());
				mData.push_back(row);
			}

			if ((columnIdx + 1) > GetDataColumnCount())
			{
				for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
				{
					itRow->resize(columnIdx + 1 + (mLabelParams.mRowNameIdx + 1));
				}
			}

			Converter<T> converter(mConverterParams);
			for (auto itRow = pColumn.begin(); itRow != pColumn.end(); ++itRow)
			{
				std::string str;
				converter.ToStr(*itRow, str);
				mData.at(std::distance(pColumn.begin(), itRow) + (mLabelParams.mColumnNameIdx + 1)).at(columnIdx) = str;
			}
		}

		/**
		 * @brief   Set column by name.
		 * @param   pColumnName           column label name.
		 * @param   pColumn               vector of column data.
		 */
		template<typename T>
		void SetColumn(const std::string& pColumnName, const std::vector<T>& pColumn)
		{
			const ssize_t columnIdx = GetColumnIdx(pColumnName);
			if (columnIdx < 0)
			{
				throw std::out_of_range("column not found: " + pColumnName);
			}
			SetColumn<T>(columnIdx, pColumn);
		}

		/**
		 * @brief   Remove column by index.
		 * @param   pColumnIdx            zero-based column index.
		 */
		void RemoveColumn(const size_t pColumnIdx)
		{
			const ssize_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
			{
				itRow->erase(itRow->begin() + columnIdx);
			}
		}

		/**
		 * @brief   Remove column by name.
		 * @param   pColumnName           column label name.
		 */
		void RemoveColumn(const std::string& pColumnName)
		{
			ssize_t columnIdx = GetColumnIdx(pColumnName);
			if (columnIdx < 0)
			{
				throw std::out_of_range("column not found: " + pColumnName);
			}

			RemoveColumn(columnIdx);
		}

		/**
		 * @brief   Get number of data columns.
		 * @returns column count.
		 */
		size_t GetColumnCount() const
		{
			return (mData.size() > 0) ? (mData.at(0).size() - (mLabelParams.mRowNameIdx + 1)) : 0;
		}

		/**
		 * @brief   Get row by index.
		 * @param   pRowIdx               zero-based row index.
		 * @returns vector of row data.
		 */
		template<typename T>
		std::vector<T> GetRow(const size_t pRowIdx) const
		{
			const ssize_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);
			std::vector<T> row;
			Converter<T> converter(mConverterParams);
			for (auto itCol = mData.at(rowIdx).begin(); itCol != mData.at(rowIdx).end(); ++itCol)
			{
				if (std::distance(mData.at(rowIdx).begin(), itCol) > mLabelParams.mRowNameIdx)
				{
					T val;
					converter.ToVal(*itCol, val);
					row.push_back(val);
				}
			}
			return row;
		}

		/**
		 * @brief   Get row by name.
		 * @param   pRowName              row label name.
		 * @returns vector of row data.
		 */
		template<typename T>
		std::vector<T> GetRow(const std::string& pRowName) const
		{
			ssize_t rowIdx = GetRowIdx(pRowName);
			if (rowIdx < 0)
			{
				throw std::out_of_range("row not found: " + pRowName);
			}
			return GetRow<T>(rowIdx);
		}

		/**
		 * @brief   Set row by index.
		 * @param   pRowIdx               zero-based row index.
		 * @param   pRow                  vector of row data.
		 */
		template<typename T>
		void SetRow(const size_t pRowIdx, const std::vector<T>& pRow)
		{
			const size_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);

			while ((rowIdx + 1) > GetDataRowCount())
			{
				std::vector<std::string> row;
				row.resize(GetDataColumnCount());
				mData.push_back(row);
			}

			if (pRow.size() > GetDataColumnCount())
			{
				for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
				{
					itRow->resize(pRow.size() + (mLabelParams.mRowNameIdx + 1));
				}
			}

			Converter<T> converter(mConverterParams);
			for (auto itCol = pRow.begin(); itCol != pRow.end(); ++itCol)
			{
				std::string str;
				converter.ToStr(*itCol, str);
				mData.at(rowIdx).at(std::distance(pRow.begin(), itCol) + (mLabelParams.mRowNameIdx + 1)) = str;
			}
		}

		/**
		 * @brief   Set row by name.
		 * @param   pRowName              row label name.
		 * @param   pRow                  vector of row data.
		 */
		template<typename T>
		void SetRow(const std::string& pRowName, const std::vector<T>& pRow)
		{
			ssize_t rowIdx = GetRowIdx(pRowName);
			if (rowIdx < 0)
			{
				throw std::out_of_range("row not found: " + pRowName);
			}
			return SetRow<T>(rowIdx, pRow);
		}

		/**
		 * @brief   Remove row by index.
		 * @param   pRowIdx               zero-based row index.
		 */
		void RemoveRow(const size_t pRowIdx)
		{
			const ssize_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);
			mData.erase(mData.begin() + rowIdx);
		}

		/**
		 * @brief   Remove row by name.
		 * @param   pRowName              row label name.
		 */
		void RemoveRow(const std::string& pRowName)
		{
			ssize_t rowIdx = GetRowIdx(pRowName);
			if (rowIdx < 0)
			{
				throw std::out_of_range("row not found: " + pRowName);
			}

			RemoveRow(rowIdx);
		}

		/**
		 * @brief   Get number of data rows.
		 * @returns row count.
		 */
		size_t GetRowCount() const
		{
			return mData.size() - (mLabelParams.mColumnNameIdx + 1);
		}

		/**
		 * @brief   Get cell by index.
		 * @param   pRowIdx               zero-based row index.
		 * @param   pColumnIdx            zero-based column index.
		 * @returns cell data.
		 */
		template<typename T>
		T GetCell(const size_t pColumnIdx, const size_t pRowIdx) const
		{
			const ssize_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			const ssize_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);

			T val;
			Converter<T> converter(mConverterParams);
			converter.ToVal(mData.at(rowIdx).at(columnIdx), val);
			return val;
		}

		/**
		 * @brief   Get cell by name.
		 * @param   pColumnName           column label name.
		 * @param   pRowName              row label name.
		 * @returns cell data.
		 */
		template<typename T>
		T GetCell(const std::string& pColumnName, const std::string& pRowName) const
		{
			const ssize_t columnIdx = GetColumnIdx(pColumnName);
			if (columnIdx < 0)
			{
				throw std::out_of_range("column not found: " + pColumnName);
			}

			const ssize_t rowIdx = GetRowIdx(pRowName);
			if (rowIdx < 0)
			{
				throw std::out_of_range("row not found: " + pRowName);
			}

			return GetCell<T>(columnIdx, rowIdx);
		}

		/**
		 * @brief   Set cell by index.
		 * @param   pRowIdx               zero-based row index.
		 * @param   pColumnIdx            zero-based column index.
		 * @param   pCell                 cell data.
		 */
		template<typename T>
		void SetCell(const size_t pColumnIdx, const size_t pRowIdx, const T& pCell)
		{
			const size_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			const size_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);

			while ((rowIdx + 1) > GetDataRowCount())
			{
				std::vector<std::string> row;
				row.resize(GetDataColumnCount());
				mData.push_back(row);
			}

			if ((columnIdx + 1) > GetDataColumnCount())
			{
				for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
				{
					itRow->resize(columnIdx + 1);
				}
			}

			std::string str;
			Converter<T> converter(mConverterParams);
			converter.ToStr(pCell, str);
			mData.at(rowIdx).at(columnIdx) = str;
		}

		/**
		 * @brief   Set cell by name.
		 * @param   pColumnName           column label name.
		 * @param   pRowName              row label name.
		 * @param   pCell                 cell data.
		 */
		template<typename T>
		void SetCell(const std::string& pColumnName, const std::string& pRowName, const T& pCell)
		{
			const ssize_t columnIdx = GetColumnIdx(pColumnName);
			if (columnIdx < 0)
			{
				throw std::out_of_range("column not found: " + pColumnName);
			}

			const ssize_t rowIdx = GetRowIdx(pRowName);
			if (rowIdx < 0)
			{
				throw std::out_of_range("row not found: " + pRowName);
			}

			SetCell<T>(columnIdx, rowIdx, pCell);
		}

		/**
		 * @brief   Get column name
		 * @param   pColumnIdx            zero-based column index.
		 * @returns column name.
		 */
		std::string GetColumnName(const ssize_t pColumnIdx)
		{
			const ssize_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			if (mLabelParams.mColumnNameIdx < 0)
			{
				throw std::out_of_range("column name row index < 0: " + std::to_string(mLabelParams.mColumnNameIdx));
			}

			return mData.at(mLabelParams.mColumnNameIdx).at(columnIdx);
		}

		/**
		 * @brief   Set column name
		 * @param   pColumnIdx            zero-based column index.
		 * @param   pColumnName           column name.
		 */
		void SetColumnName(size_t pColumnIdx, const std::string& pColumnName)
		{
			const ssize_t columnIdx = pColumnIdx + (mLabelParams.mRowNameIdx + 1);
			mColumnNames[pColumnName] = columnIdx;
			if (mLabelParams.mColumnNameIdx < 0)
			{
				throw std::out_of_range("column name row index < 0: " + std::to_string(mLabelParams.mColumnNameIdx));
			}

			mData.at(mLabelParams.mColumnNameIdx).at(columnIdx) = pColumnName;
		}

		/**
		 * @brief   Get column names
		 * @returns vector of column names.
		 */
		std::vector<std::string> GetColumnNames()
		{
			if (mLabelParams.mColumnNameIdx >= 0)
			{
				return std::vector<std::string>(mData.at(mLabelParams.mColumnNameIdx).begin() +
					(mLabelParams.mRowNameIdx + 1),
					mData.at(mLabelParams.mColumnNameIdx).end());
			}

			return std::vector<std::string>();
		}

		/**
		 * @brief   Get row name
		 * @param   pRowIdx               zero-based column index.
		 * @returns row name.
		 */
		std::string GetRowName(const ssize_t pRowIdx)
		{
			const ssize_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);
			if (mLabelParams.mRowNameIdx < 0)
			{
				throw std::out_of_range("row name column index < 0: " + std::to_string(mLabelParams.mRowNameIdx));
			}

			return mData.at(rowIdx).at(mLabelParams.mRowNameIdx);
		}

		/**
		 * @brief   Set row name
		 * @param   pRowIdx               zero-based row index.
		 * @param   pRowName              row name.
		 */
		void SetRowName(size_t pRowIdx, const std::string& pRowName)
		{
			const ssize_t rowIdx = pRowIdx + (mLabelParams.mColumnNameIdx + 1);
			mRowNames[pRowName] = rowIdx;
			if (mLabelParams.mRowNameIdx < 0)
			{
				throw std::out_of_range("row name column index < 0: " + std::to_string(mLabelParams.mRowNameIdx));
			}

			mData.at(rowIdx).at(mLabelParams.mRowNameIdx) = pRowName;
		}

		/**
		 * @brief   Get row names
		 * @returns vector of row names.
		 */
		std::vector<std::string> GetRowNames()
		{
			std::vector<std::string> rownames;
			if (mLabelParams.mRowNameIdx >= 0)
			{
				for (auto itRow = mData.begin(); itRow != mData.end(); ++itRow)
				{
					if (std::distance(mData.begin(), itRow) > mLabelParams.mColumnNameIdx)
					{
						rownames.push_back(itRow->at(mLabelParams.mRowNameIdx));
					}
				}
			}
			return rownames;
		}

	private:
		void ReadCsv()
		{
			std::ifstream stream;
			stream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
			stream.open(mPath, std::ios::binary);

#ifdef HAS_CODECVT
			stream.seekg(0, std::ios::end);
			std::streamsize length = stream.tellg();
			stream.seekg(0, std::ios::beg);

			std::vector<char> bom(2, '\0');
			if (length >= 2)
			{
				stream.read(bom.data(), 2);
			}

			static const std::vector<char> bomU16le = { '\xff', '\xfe' };
			static const std::vector<char> bomU16be = { '\xfe', '\xff' };
			if ((bom == bomU16le) || (bom == bomU16be))
			{
				mIsUtf16 = true;
				mIsLE = (bom == bomU16le);

				std::wifstream wstream;
				wstream.exceptions(std::wifstream::failbit | std::wifstream::badbit);
				wstream.open(mPath, std::ios::binary);
				if (mIsLE)
				{
					wstream.imbue(std::locale(wstream.getloc(),
						new std::codecvt_utf16<wchar_t, 0x10ffff,
						static_cast<std::codecvt_mode>(std::consume_header | std::little_endian)>));
				}
				else
				{
					wstream.imbue(std::locale(wstream.getloc(),
						new std::codecvt_utf16<wchar_t, 0x10ffff,
						std::consume_header>));
				}
				std::wstringstream wss;
				wss << wstream.rdbuf();
				std::string utf8 = ToString(wss.str());
				std::stringstream ss(utf8);
				ReadCsv(ss);
			}
			else
#endif
			{
				stream.seekg(0, std::ios::beg);
				ReadCsv(stream);
			}
		}

		void ReadCsv(std::istream& pStream)
		{
			pStream.seekg(0, std::ios::end);
			std::streamsize fileLength = pStream.tellg();
			pStream.seekg(0, std::ios::beg);
			const std::streamsize bufLength = 64 * 1024;
			std::vector<char> buffer(bufLength);
			std::vector<std::string> row;
			std::string cell;
			bool quoted = false;
			int cr = 0;
			int lf = 0;

			while (fileLength > 0)
			{
				std::streamsize readLength = std::min(fileLength, bufLength);
				pStream.read(buffer.data(), readLength);
				for (int i = 0; i < readLength; ++i)
				{
					if (buffer[i] == '"')
					{
						if (cell.empty() || cell[0] == '"')
						{
							quoted = !quoted;
						}
						cell += buffer[i];
					}
					else if (buffer[i] == mSeparatorParams.mSeparator)
					{
						if (!quoted)
						{
							row.push_back(cell);
							cell.clear();
						}
						else
						{
							cell += buffer[i];
						}
					}
					else if (buffer[i] == '\r')
					{
						++cr;
					}
					else if (buffer[i] == '\n')
					{
						++lf;
						row.push_back(cell);
						cell.clear();
						mData.push_back(row);
						row.clear();
						quoted = false; // disallow line breaks in quoted string, by auto-unquote at linebreak
					}
					else
					{
						cell += buffer[i];
					}
				}
				fileLength -= readLength;
			}

			// Handle last line without linebreak
			if (!cell.empty() || !row.empty())
			{
				row.push_back(cell);
				cell.clear();
				mData.push_back(row);
				row.clear();
			}

			// Assume CR/LF if at least half the linebreaks have CR
			mSeparatorParams.mHasCR = (cr > (lf / 2));

			// Set up column labels
			if ((mLabelParams.mColumnNameIdx >= 0) &&
				(mData.size() > 0))
			{
				int i = 0;
				for (auto& columnName : mData[mLabelParams.mColumnNameIdx])
				{
					mColumnNames[columnName] = i++;
				}
			}

			// Set up row labels
			if ((mLabelParams.mRowNameIdx >= 0) &&
				(static_cast<ssize_t>(mData.size()) >
				(mLabelParams.mColumnNameIdx + 1)))
			{
				int i = 0;
				for (auto& dataRow : mData)
				{
					mRowNames[dataRow[mLabelParams.mRowNameIdx]] = i++;
				}
			}
		}

		void WriteCsv() const
		{
#ifdef HAS_CODECVT
			if (mIsUtf16)
			{
				std::stringstream ss;
				WriteCsv(ss);
				std::string utf8 = ss.str();
				std::wstring wstr = ToWString(utf8);

				std::wofstream wstream;
				wstream.exceptions(std::wofstream::failbit | std::wofstream::badbit);
				wstream.open(mPath, std::ios::binary | std::ios::trunc);

				if (mIsLE)
				{
					wstream.imbue(std::locale(wstream.getloc(),
						new std::codecvt_utf16<wchar_t, 0x10ffff,
						static_cast<std::codecvt_mode>(std::little_endian)>));
				}
				else
				{
					wstream.imbue(std::locale(wstream.getloc(),
						new std::codecvt_utf16<wchar_t, 0x10ffff>));
				}

				wstream << (wchar_t)0xfeff;
				wstream << wstr;
			}
			else
#endif
			{
				std::ofstream stream;
				stream.exceptions(std::ofstream::failbit | std::ofstream::badbit);
				stream.open(mPath, std::ios::binary | std::ios::trunc);
				WriteCsv(stream);
			}
		}

		void WriteCsv(std::ostream& pStream) const
		{
			for (auto itr = mData.begin(); itr != mData.end(); ++itr)
			{
				for (auto itc = itr->begin(); itc != itr->end(); ++itc)
				{
					if ((std::string::npos == itc->find(mSeparatorParams.mSeparator)) ||
						((itc->length() >= 2) && ((*itc)[0] == '\"') && ((*itc)[itc->length() - 1] == '\"')))
					{
						pStream << *itc;
					}
					else
					{
						pStream << '"' << *itc << '"';
					}

					if (std::distance(itc, itr->end()) > 1)
					{
						pStream << mSeparatorParams.mSeparator;
					}
				}
				pStream << (mSeparatorParams.mHasCR ? "\r\n" : "\n");
			}
		}

		ssize_t GetColumnIdx(const std::string& pColumnName) const
		{
			if (mLabelParams.mColumnNameIdx >= 0)
			{
				if (mColumnNames.find(pColumnName) != mColumnNames.end())
				{
					return mColumnNames.at(pColumnName) - (mLabelParams.mRowNameIdx + 1);
				}
			}
			return -1;
		}

		ssize_t GetRowIdx(const std::string& pRowName) const
		{
			if (mLabelParams.mRowNameIdx >= 0)
			{
				if (mRowNames.find(pRowName) != mRowNames.end())
				{
					return mRowNames.at(pRowName) - (mLabelParams.mColumnNameIdx + 1);
				}
			}
			return -1;
		}

		size_t GetDataRowCount() const
		{
			return mData.size();
		}

		size_t GetDataColumnCount() const
		{
			return (mData.size() > 0) ? mData.at(0).size() : 0;
		}

#ifdef HAS_CODECVT
#if defined(_MSC_VER)
#pragma warning (disable: 4996)
#endif
		static std::string ToString(const std::wstring& pWStr)
		{
			size_t len = std::wcstombs(nullptr, pWStr.c_str(), 0) + 1;
			char* cstr = new char[len];
			std::wcstombs(cstr, pWStr.c_str(), len);
			std::string str(cstr);
			delete[] cstr;
			return str;
		}

		static std::wstring ToWString(const std::string& pStr)
		{
			size_t len = 1 + mbstowcs(nullptr, pStr.c_str(), 0);
			wchar_t* wcstr = new wchar_t[len];
			std::mbstowcs(wcstr, pStr.c_str(), len);
			std::wstring wstr(wcstr);
			delete[] wcstr;
			return wstr;
		}
#if defined(_MSC_VER)
#pragma warning (default: 4996)
#endif
#endif

	private:
		std::string mPath;
		LabelParams mLabelParams;
		SeparatorParams mSeparatorParams;
		ConverterParams mConverterParams;
		std::vector<std::vector<std::string> > mData;
		std::map<std::string, size_t> mColumnNames;
		std::map<std::string, size_t> mRowNames;
#ifdef HAS_CODECVT
		bool mIsUtf16 = false;
		bool mIsLE = false;
#endif
	};
}


namespace io {
	////////////////////////////////////////////////////////////////////////////
	//                                 LineReader                             //
	////////////////////////////////////////////////////////////////////////////

	namespace error {
		struct base : std::exception {
			virtual void format_error_message()const = 0;

			const char*what()const throw() {
				format_error_message();
				return error_message_buffer;
			}

			mutable char error_message_buffer[512];
		};

		const int max_file_name_length = 255;

		struct with_file_name {
			with_file_name() {
				std::memset(file_name, 0, sizeof(file_name));
			}

			void set_file_name(const char*file_name) {
				if (file_name != nullptr) {
					strncpy(this->file_name, file_name, sizeof(this->file_name));
					this->file_name[sizeof(this->file_name) - 1] = '\0';
				}
				else {
					this->file_name[0] = '\0';
				}
			}

			char file_name[max_file_name_length + 1];
		};

		struct with_file_line {
			with_file_line() {
				file_line = -1;
			}

			void set_file_line(int file_line) {
				this->file_line = file_line;
			}

			int file_line;
		};

		struct with_errno {
			with_errno() {
				errno_value = 0;
			}

			void set_errno(int errno_value) {
				this->errno_value = errno_value;
			}

			int errno_value;
		};

		struct can_not_open_file :
			base,
			with_file_name,
			with_errno {
			void format_error_message()const {
				if (errno_value != 0)
					std::snprintf(error_message_buffer, sizeof(error_message_buffer),
						"Can not open file \"%s\" because \"%s\"."
						, file_name, std::strerror(errno_value));
				else
					std::snprintf(error_message_buffer, sizeof(error_message_buffer),
						"Can not open file \"%s\"."
						, file_name);
			}
		};

		struct line_length_limit_exceeded :
			base,
			with_file_name,
			with_file_line {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Line number %d in file \"%s\" exceeds the maximum length of 2^24-1."
					, file_line, file_name);
			}
		};
	}

	class ByteSourceBase {
	public:
		virtual int read(char*buffer, int size) = 0;
		virtual ~ByteSourceBase() {}
	};

	namespace detail {
		class OwningStdIOByteSourceBase : public ByteSourceBase {
		public:
			explicit OwningStdIOByteSourceBase(FILE*file) :file(file) {
				// Tell the std library that we want to do the buffering ourself.
				std::setvbuf(file, 0, _IONBF, 0);
			}

			int read(char*buffer, int size) {
				return std::fread(buffer, 1, size, file);
			}

			~OwningStdIOByteSourceBase() {
				std::fclose(file);
			}

		private:
			FILE*file;
		};

		class NonOwningIStreamByteSource : public ByteSourceBase {
		public:
			explicit NonOwningIStreamByteSource(std::istream&in) :in(in) {}

			int read(char*buffer, int size) {
				in.read(buffer, size);
				return in.gcount();
			}

			~NonOwningIStreamByteSource() {}

		private:
			std::istream&in;
		};

		class NonOwningStringByteSource : public ByteSourceBase {
		public:
			NonOwningStringByteSource(const char*str, long long size) :str(str), remaining_byte_count(size) {}

			int read(char*buffer, int desired_byte_count) {
				int to_copy_byte_count = desired_byte_count;
				if (remaining_byte_count < to_copy_byte_count)
					to_copy_byte_count = remaining_byte_count;
				std::memcpy(buffer, str, to_copy_byte_count);
				remaining_byte_count -= to_copy_byte_count;
				str += to_copy_byte_count;
				return to_copy_byte_count;
			}

			~NonOwningStringByteSource() {}

		private:
			const char*str;
			long long remaining_byte_count;
		};

		class SynchronousReader {
		public:
			void init(std::unique_ptr<ByteSourceBase>arg_byte_source) {
				byte_source = std::move(arg_byte_source);
			}

			bool is_valid()const {
				return byte_source != nullptr;
			}

			void start_read(char*arg_buffer, int arg_desired_byte_count) {
				buffer = arg_buffer;
				desired_byte_count = arg_desired_byte_count;
			}

			int finish_read() {
				return byte_source->read(buffer, desired_byte_count);
			}
		private:
			std::unique_ptr<ByteSourceBase>byte_source;
			char*buffer;
			int desired_byte_count;
		};
	}

	class LineReader {
	private:
		static const int block_len = 1 << 24;
		std::unique_ptr<char[]>buffer; // must be constructed before (and thus destructed after) the reader!
		detail::SynchronousReader reader;
		int data_begin;
		int data_end;

		char file_name[error::max_file_name_length + 1];
		unsigned file_line;

		static std::unique_ptr<ByteSourceBase> open_file(const char*file_name) {
			// We open the file in binary mode as it makes no difference under *nix
			// and under Windows we handle \r\n newlines ourself.
			FILE*file = std::fopen(file_name, "rb");
			if (file == 0) {
				int x = errno; // store errno as soon as possible, doing it after constructor call can fail.
				error::can_not_open_file err;
				err.set_errno(x);
				err.set_file_name(file_name);
				throw err;
			}
			return std::unique_ptr<ByteSourceBase>(new detail::OwningStdIOByteSourceBase(file));
		}

		void init(std::unique_ptr<ByteSourceBase>byte_source) {
			file_line = 0;

			buffer = std::unique_ptr<char[]>(new char[3 * block_len]);
			data_begin = 0;
			data_end = byte_source->read(buffer.get(), 2 * block_len);

			// Ignore UTF-8 BOM
			if (data_end >= 3 && buffer[0] == '\xEF' && buffer[1] == '\xBB' && buffer[2] == '\xBF')
				data_begin = 3;

			if (data_end == 2 * block_len) {
				reader.init(std::move(byte_source));
				reader.start_read(buffer.get() + 2 * block_len, block_len);
			}
		}

	public:
		LineReader() = delete;
		LineReader(const LineReader&) = delete;
		LineReader&operator=(const LineReader&) = delete;

		explicit LineReader(const char*file_name) {
			set_file_name(file_name);
			init(open_file(file_name));
		}

		explicit LineReader(const std::string&file_name) {
			set_file_name(file_name.c_str());
			init(open_file(file_name.c_str()));
		}

		LineReader(const char*file_name, std::unique_ptr<ByteSourceBase>byte_source) {
			set_file_name(file_name);
			init(std::move(byte_source));
		}

		LineReader(const std::string&file_name, std::unique_ptr<ByteSourceBase>byte_source) {
			set_file_name(file_name.c_str());
			init(std::move(byte_source));
		}

		LineReader(const char*file_name, const char*data_begin, const char*data_end) {
			set_file_name(file_name);
			init(std::unique_ptr<ByteSourceBase>(new detail::NonOwningStringByteSource(data_begin, data_end - data_begin)));
		}

		LineReader(const std::string&file_name, const char*data_begin, const char*data_end) {
			set_file_name(file_name.c_str());
			init(std::unique_ptr<ByteSourceBase>(new detail::NonOwningStringByteSource(data_begin, data_end - data_begin)));
		}

		LineReader(const char*file_name, FILE*file) {
			set_file_name(file_name);
			init(std::unique_ptr<ByteSourceBase>(new detail::OwningStdIOByteSourceBase(file)));
		}

		LineReader(const std::string&file_name, FILE*file) {
			set_file_name(file_name.c_str());
			init(std::unique_ptr<ByteSourceBase>(new detail::OwningStdIOByteSourceBase(file)));
		}

		LineReader(const char*file_name, std::istream&in) {
			set_file_name(file_name);
			init(std::unique_ptr<ByteSourceBase>(new detail::NonOwningIStreamByteSource(in)));
		}

		LineReader(const std::string&file_name, std::istream&in) {
			set_file_name(file_name.c_str());
			init(std::unique_ptr<ByteSourceBase>(new detail::NonOwningIStreamByteSource(in)));
		}

		void set_file_name(const std::string&file_name) {
			set_file_name(file_name.c_str());
		}

		void set_file_name(const char*file_name) {
			if (file_name != nullptr) {
				strncpy(this->file_name, file_name, sizeof(this->file_name));
				this->file_name[sizeof(this->file_name) - 1] = '\0';
			}
			else {
				this->file_name[0] = '\0';
			}
		}

		const char*get_truncated_file_name()const {
			return file_name;
		}

		void set_file_line(unsigned file_line) {
			this->file_line = file_line;
		}

		unsigned get_file_line()const {
			return file_line;
		}

		char*next_line() {
			if (data_begin == data_end)
				return 0;

			++file_line;

			assert(data_begin < data_end);
			assert(data_end <= block_len * 2);

			if (data_begin >= block_len) {
				std::memcpy(buffer.get(), buffer.get() + block_len, block_len);
				data_begin -= block_len;
				data_end -= block_len;
				if (reader.is_valid())
				{
					data_end += reader.finish_read();
					std::memcpy(buffer.get() + block_len, buffer.get() + 2 * block_len, block_len);
					reader.start_read(buffer.get() + 2 * block_len, block_len);
				}
			}

			int line_end = data_begin;
			while (buffer[line_end] != '\n' && line_end != data_end) {
				++line_end;
			}

			if (line_end - data_begin + 1 > block_len) {
				error::line_length_limit_exceeded err;
				err.set_file_name(file_name);
				err.set_file_line(file_line);
				throw err;
			}

			if (buffer[line_end] == '\n' && line_end != data_end) {
				buffer[line_end] = '\0';
			}
			else {
				// some files are missing the newline at the end of the
				// last line
				++data_end;
				buffer[line_end] = '\0';
			}

			// handle windows \r\n-line breaks
			if (line_end != data_begin && buffer[line_end - 1] == '\r')
				buffer[line_end - 1] = '\0';

			char*ret = buffer.get() + data_begin;
			data_begin = line_end + 1;
			return ret;
		}
	};

	////////////////////////////////////////////////////////////////////////////
	//                                 CSV                                    //
	////////////////////////////////////////////////////////////////////////////

	namespace error {
		const int max_column_name_length = 63;
		struct with_column_name {
			with_column_name() {
				std::memset(column_name, 0, max_column_name_length + 1);
			}

			void set_column_name(const char*column_name) {
				if (column_name != nullptr) {
					std::strncpy(this->column_name, column_name, max_column_name_length);
					this->column_name[max_column_name_length] = '\0';
				}
				else {
					this->column_name[0] = '\0';
				}
			}

			char column_name[max_column_name_length + 1];
		};

		const int max_column_content_length = 63;

		struct with_column_content {
			with_column_content() {
				std::memset(column_content, 0, max_column_content_length + 1);
			}

			void set_column_content(const char*column_content) {
				if (column_content != nullptr) {
					std::strncpy(this->column_content, column_content, max_column_content_length);
					this->column_content[max_column_content_length] = '\0';
				}
				else {
					this->column_content[0] = '\0';
				}
			}

			char column_content[max_column_content_length + 1];
		};

		struct extra_column_in_header :
			base,
			with_file_name,
			with_column_name {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Extra column \"%s\" in header of file \"%s\"."
					, column_name, file_name);
			}
		};

		struct missing_column_in_header :
			base,
			with_file_name,
			with_column_name {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Missing column \"%s\" in header of file \"%s\"."
					, column_name, file_name);
			}
		};

		struct duplicated_column_in_header :
			base,
			with_file_name,
			with_column_name {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Duplicated column \"%s\" in header of file \"%s\"."
					, column_name, file_name);
			}
		};

		struct header_missing :
			base,
			with_file_name {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Header missing in file \"%s\"."
					, file_name);
			}
		};

		struct too_few_columns :
			base,
			with_file_name,
			with_file_line {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Too few columns in line %d in file \"%s\"."
					, file_line, file_name);
			}
		};

		struct too_many_columns :
			base,
			with_file_name,
			with_file_line {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Too many columns in line %d in file \"%s\"."
					, file_line, file_name);
			}
		};

		struct escaped_string_not_closed :
			base,
			with_file_name,
			with_file_line {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"Escaped string was not closed in line %d in file \"%s\"."
					, file_line, file_name);
			}
		};

		struct integer_must_be_positive :
			base,
			with_file_name,
			with_file_line,
			with_column_name,
			with_column_content {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"The integer \"%s\" must be positive or 0 in column \"%s\" in file \"%s\" in line \"%d\"."
					, column_content, column_name, file_name, file_line);
			}
		};

		struct no_digit :
			base,
			with_file_name,
			with_file_line,
			with_column_name,
			with_column_content {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"The integer \"%s\" contains an invalid digit in column \"%s\" in file \"%s\" in line \"%d\"."
					, column_content, column_name, file_name, file_line);
			}
		};

		struct integer_overflow :
			base,
			with_file_name,
			with_file_line,
			with_column_name,
			with_column_content {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"The integer \"%s\" overflows in column \"%s\" in file \"%s\" in line \"%d\"."
					, column_content, column_name, file_name, file_line);
			}
		};

		struct integer_underflow :
			base,
			with_file_name,
			with_file_line,
			with_column_name,
			with_column_content {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"The integer \"%s\" underflows in column \"%s\" in file \"%s\" in line \"%d\"."
					, column_content, column_name, file_name, file_line);
			}
		};

		struct invalid_single_character :
			base,
			with_file_name,
			with_file_line,
			with_column_name,
			with_column_content {
			void format_error_message()const {
				std::snprintf(error_message_buffer, sizeof(error_message_buffer),
					"The content \"%s\" of column \"%s\" in file \"%s\" in line \"%d\" is not a single character."
					, column_content, column_name, file_name, file_line);
			}
		};
	}

	typedef unsigned ignore_column;
	static const ignore_column ignore_no_column = 0;
	static const ignore_column ignore_extra_column = 1;
	static const ignore_column ignore_missing_column = 2;

	template<char ... trim_char_list>
	struct trim_chars {
	private:
		constexpr static bool is_trim_char(char) {
			return false;
		}

		template<class ...OtherTrimChars>
		constexpr static bool is_trim_char(char c, char trim_char, OtherTrimChars...other_trim_chars) {
			return c == trim_char || is_trim_char(c, other_trim_chars...);
		}

	public:
		static void trim(char*&str_begin, char*&str_end) {
			while (str_begin != str_end && is_trim_char(*str_begin, trim_char_list...))
				++str_begin;
			while (str_begin != str_end && is_trim_char(*(str_end - 1), trim_char_list...))
				--str_end;
			*str_end = '\0';
		}
	};

	struct no_comment {
		static bool is_comment(const char*) {
			return false;
		}
	};

	template<char ... comment_start_char_list>
	struct single_line_comment {
	private:
		constexpr static bool is_comment_start_char(char) {
			return false;
		}

		template<class ...OtherCommentStartChars>
		constexpr static bool is_comment_start_char(char c, char comment_start_char, OtherCommentStartChars...other_comment_start_chars) {
			return c == comment_start_char || is_comment_start_char(c, other_comment_start_chars...);
		}

	public:

		static bool is_comment(const char*line) {
			return is_comment_start_char(*line, comment_start_char_list...);
		}
	};

	struct empty_line_comment {
		static bool is_comment(const char*line) {
			if (*line == '\0')
				return true;
			while (*line == ' ' || *line == '\t') {
				++line;
				if (*line == 0)
					return true;
			}
			return false;
		}
	};

	template<char ... comment_start_char_list>
	struct single_and_empty_line_comment {
		static bool is_comment(const char*line) {
			return single_line_comment<comment_start_char_list...>::is_comment(line) || empty_line_comment::is_comment(line);
		}
	};

	template<char sep>
	struct no_quote_escape {
		static const char*find_next_column_end(const char*col_begin) {
			while (*col_begin != sep && *col_begin != '\0')
				++col_begin;
			return col_begin;
		}

		static void unescape(char*&, char*&) {
		}
	};

	template<char sep, char quote>
	struct double_quote_escape {
		static const char*find_next_column_end(const char*col_begin) {
			while (*col_begin != sep && *col_begin != '\0')
				if (*col_begin != quote)
					++col_begin;
				else {
					do {
						++col_begin;
						while (*col_begin != quote) {
							if (*col_begin == '\0')
								throw error::escaped_string_not_closed();
							++col_begin;
						}
						++col_begin;
					} while (*col_begin == quote);
				}
			return col_begin;
		}

		static void unescape(char*&col_begin, char*&col_end) {
			if (col_end - col_begin >= 2) {
				if (*col_begin == quote && *(col_end - 1) == quote) {
					++col_begin;
					--col_end;
					char*out = col_begin;
					for (char*in = col_begin; in != col_end; ++in) {
						if (*in == quote && (in + 1) != col_end && *(in + 1) == quote) {
							++in;
						}
						*out = *in;
						++out;
					}
					col_end = out;
					*col_end = '\0';
				}
			}
		}
	};

	struct throw_on_overflow {
		template<class T>
		static void on_overflow(T&) {
			throw error::integer_overflow();
		}

		template<class T>
		static void on_underflow(T&) {
			throw error::integer_underflow();
		}
	};

	struct ignore_overflow {
		template<class T>
		static void on_overflow(T&) {}

		template<class T>
		static void on_underflow(T&) {}
	};

	struct set_to_max_on_overflow {
		template<class T>
		static void on_overflow(T&x) {
			x = std::numeric_limits<T>::max();
		}

		template<class T>
		static void on_underflow(T&x) {
			x = std::numeric_limits<T>::min();
		}
	};

	namespace detail {
		template<class quote_policy>
		void chop_next_column(
			char*&line, char*&col_begin, char*&col_end
		) {
			assert(line != nullptr);

			col_begin = line;
			// the col_begin + (... - col_begin) removes the constness
			col_end = col_begin + (quote_policy::find_next_column_end(col_begin) - col_begin);

			if (*col_end == '\0') {
				line = nullptr;
			}
			else {
				*col_end = '\0';
				line = col_end + 1;
			}
		}

		template<class trim_policy, class quote_policy>
		void parse_line(
			char*line,
			char**sorted_col,
			const std::vector<int>&col_order
		) {
			for (std::size_t i = 0; i < col_order.size(); ++i) {
				if (line == nullptr)
					throw ::io::error::too_few_columns();
				char*col_begin, *col_end;
				chop_next_column<quote_policy>(line, col_begin, col_end);

				if (col_order[i] != -1) {
					trim_policy::trim(col_begin, col_end);
					quote_policy::unescape(col_begin, col_end);

					sorted_col[col_order[i]] = col_begin;
				}
			}
			if (line != nullptr)
				throw ::io::error::too_many_columns();
		}

		template<unsigned column_count, class trim_policy, class quote_policy>
		void parse_header_line(
			char*line,
			std::vector<int>&col_order,
			const std::string*col_name,
			ignore_column ignore_policy
		) {
			col_order.clear();

			bool found[column_count];
			std::fill(found, found + column_count, false);
			while (line) {
				char*col_begin, *col_end;
				chop_next_column<quote_policy>(line, col_begin, col_end);

				trim_policy::trim(col_begin, col_end);
				quote_policy::unescape(col_begin, col_end);

				for (unsigned i = 0; i < column_count; ++i)
					if (col_begin == col_name[i]) {
						if (found[i]) {
							error::duplicated_column_in_header err;
							err.set_column_name(col_begin);
							throw err;
						}
						found[i] = true;
						col_order.push_back(i);
						col_begin = 0;
						break;
					}
				if (col_begin) {
					if (ignore_policy & ::io::ignore_extra_column)
						col_order.push_back(-1);
					else {
						error::extra_column_in_header err;
						err.set_column_name(col_begin);
						throw err;
					}
				}
			}
			if (!(ignore_policy & ::io::ignore_missing_column)) {
				for (unsigned i = 0; i < column_count; ++i) {
					if (!found[i]) {
						error::missing_column_in_header err;
						err.set_column_name(col_name[i].c_str());
						throw err;
					}
				}
			}
		}

		template<class overflow_policy>
		void parse(char*col, char &x) {
			if (!*col)
				throw error::invalid_single_character();
			x = *col;
			++col;
			if (*col)
				throw error::invalid_single_character();
		}

		template<class overflow_policy>
		void parse(char*col, std::string&x) {
			x = col;
		}

		template<class overflow_policy>
		void parse(char*col, const char*&x) {
			x = col;
		}

		template<class overflow_policy>
		void parse(char*col, char*&x) {
			x = col;
		}

		template<class overflow_policy, class T>
		void parse_unsigned_integer(const char*col, T&x) {
			x = 0;
			while (*col != '\0') {
				if ('0' <= *col && *col <= '9') {
					T y = *col - '0';
					if (x > (std::numeric_limits<T>::max() - y) / 10) {
						overflow_policy::on_overflow(x);
						return;
					}
					x = 10 * x + y;
				}
				else
					throw error::no_digit();
				++col;
			}
		}

		template<class overflow_policy>void parse(char*col, unsigned char &x)
		{
			parse_unsigned_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, unsigned short &x)
		{
			parse_unsigned_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, unsigned int &x)
		{
			parse_unsigned_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, unsigned long &x)
		{
			parse_unsigned_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, unsigned long long &x)
		{
			parse_unsigned_integer<overflow_policy>(col, x);
		}

		template<class overflow_policy, class T>
		void parse_signed_integer(const char*col, T&x) {
			if (*col == '-') {
				++col;

				x = 0;
				while (*col != '\0') {
					if ('0' <= *col && *col <= '9') {
						T y = *col - '0';
						if (x < (std::numeric_limits<T>::min() + y) / 10) {
							overflow_policy::on_underflow(x);
							return;
						}
						x = 10 * x - y;
					}
					else
						throw error::no_digit();
					++col;
				}
				return;
			}
			else if (*col == '+')
				++col;
			parse_unsigned_integer<overflow_policy>(col, x);
		}

		template<class overflow_policy>void parse(char*col, signed char &x)
		{
			parse_signed_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, signed short &x)
		{
			parse_signed_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, signed int &x)
		{
			parse_signed_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, signed long &x)
		{
			parse_signed_integer<overflow_policy>(col, x);
		}
		template<class overflow_policy>void parse(char*col, signed long long &x)
		{
			parse_signed_integer<overflow_policy>(col, x);
		}

		template<class T>
		void parse_float(const char*col, T&x) {
			bool is_neg = false;
			if (*col == '-') {
				is_neg = true;
				++col;
			}
			else if (*col == '+')
				++col;

			x = 0;
			while ('0' <= *col && *col <= '9') {
				int y = *col - '0';
				x *= 10;
				x += y;
				++col;
			}

			if (*col == '.' || *col == ',') {
				++col;
				T pos = 1;
				while ('0' <= *col && *col <= '9') {
					pos /= 10;
					int y = *col - '0';
					++col;
					x += y * pos;
				}
			}

			if (*col == 'e' || *col == 'E') {
				++col;
				int e;

				parse_signed_integer<set_to_max_on_overflow>(col, e);

				if (e != 0) {
					T base;
					if (e < 0) {
						base = 0.1;
						e = -e;
					}
					else {
						base = 10;
					}

					while (e != 1) {
						if ((e & 1) == 0) {
							base = base * base;
							e >>= 1;
						}
						else {
							x *= base;
							--e;
						}
					}
					x *= base;
				}
			}
			else {
				if (*col != '\0')
					throw error::no_digit();
			}

			if (is_neg)
				x = -x;
		}

		template<class overflow_policy> void parse(char*col, float&x) { parse_float(col, x); }
		template<class overflow_policy> void parse(char*col, double&x) { parse_float(col, x); }
		template<class overflow_policy> void parse(char*col, long double&x) { parse_float(col, x); }

		template<class overflow_policy, class T>
		void parse(char*col, T&x) {
			// Mute unused variable compiler warning
			(void)col;
			(void)x;
			// GCC evalutes "false" when reading the template and
			// "sizeof(T)!=sizeof(T)" only when instantiating it. This is why
			// this strange construct is used.
			static_assert(sizeof(T) != sizeof(T),
				"Can not parse this type. Only buildin integrals, floats, char, char*, const char* and std::string are supported");
		}
	}

	template<unsigned column_count,
		class trim_policy = trim_chars<' ', '\t'>,
		class quote_policy = no_quote_escape<','>,
		class overflow_policy = throw_on_overflow,
		class comment_policy = no_comment
	>
		class CSVReader {
		private:
			LineReader in;

			char*row[column_count];
			std::string column_names[column_count];

			std::vector<int>col_order;

			template<class ...ColNames>
			void set_column_names(std::string s, ColNames...cols) {
				column_names[column_count - sizeof...(ColNames) - 1] = std::move(s);
				set_column_names(std::forward<ColNames>(cols)...);
			}

			void set_column_names() {}

		public:
			CSVReader() = delete;
			CSVReader(const CSVReader&) = delete;
			CSVReader&operator=(const CSVReader&);

			template<class ...Args>
			explicit CSVReader(Args&&...args) :in(std::forward<Args>(args)...) {
				std::fill(row, row + column_count, nullptr);
				col_order.resize(column_count);
				for (unsigned i = 0; i < column_count; ++i)
					col_order[i] = i;
				for (unsigned i = 1; i <= column_count; ++i)
					column_names[i - 1] = "col" + std::to_string(i);
			}

			char*next_line() {
				return in.next_line();
			}

			template<class ...ColNames>
			void read_header(ignore_column ignore_policy, ColNames...cols) {
				// 				static_assert(sizeof...(ColNames) >= column_count, "not enough column names specified");
				// 				static_assert(sizeof...(ColNames) <= column_count, "too many column names specified");
				try {
					set_column_names(std::forward<ColNames>(cols)...);

					char*line;
					do {
						line = in.next_line();
						if (!line)
							throw error::header_missing();
					} while (comment_policy::is_comment(line));

					detail::parse_header_line
						<column_count, trim_policy, quote_policy>
						(line, col_order, column_names, ignore_policy);
				}
				catch (error::with_file_name&err) {
					err.set_file_name(in.get_truncated_file_name());
					throw;
				}
			}

			template<class ...ColNames>
			void set_header(ColNames...cols) {
				static_assert(sizeof...(ColNames) >= column_count,
					"not enough column names specified");
				static_assert(sizeof...(ColNames) <= column_count,
					"too many column names specified");
				set_column_names(std::forward<ColNames>(cols)...);
				std::fill(row, row + column_count, nullptr);
				col_order.resize(column_count);
				for (unsigned i = 0; i < column_count; ++i)
					col_order[i] = i;
			}

			bool has_column(const std::string&name) const {
				return col_order.end() != std::find(
					col_order.begin(), col_order.end(),
					std::find(std::begin(column_names), std::end(column_names), name)
					- std::begin(column_names));
			}

			void set_file_name(const std::string&file_name) {
				in.set_file_name(file_name);
			}

			void set_file_name(const char*file_name) {
				in.set_file_name(file_name);
			}

			const char*get_truncated_file_name()const {
				return in.get_truncated_file_name();
			}

			void set_file_line(unsigned file_line) {
				in.set_file_line(file_line);
			}

			unsigned get_file_line()const {
				return in.get_file_line();
			}

		private:
			void parse_helper(std::size_t) {}

			template<class T, class ...ColType>
			void parse_helper(std::size_t r, T&t, ColType&...cols) {
				if (row[r]) {
					try {
						try {
							::io::detail::parse<overflow_policy>(row[r], t);
						}
						catch (error::with_column_content&err) {
							err.set_column_content(row[r]);
							throw;
						}
					}
					catch (error::with_column_name&err) {
						err.set_column_name(column_names[r].c_str());
						throw;
					}
				}
				parse_helper(r + 1, cols...);
			}

		public:
			template<class ...ColType>
			bool read_row(ColType& ...cols) {
				// 				static_assert(sizeof...(ColType) >= column_count,
				// 					"not enough columns specified");
				// 				static_assert(sizeof...(ColType) <= column_count,
				// 					"too many columns specified");
				try {
					try {
						char*line;
						do {
							line = in.next_line();
							if (!line)
								return false;
						} while (comment_policy::is_comment(line));

						detail::parse_line<trim_policy, quote_policy>
							(line, row, col_order);

						parse_helper(0, cols...);
					}
					catch (error::with_file_name&err) {
						err.set_file_name(in.get_truncated_file_name());
						throw;
					}
				}
				catch (error::with_file_line&err) {
					err.set_file_line(in.get_file_line());
					throw;
				}

				return true;
			}
	};
}
#endif