import { motion } from 'framer-motion'
import { ArrowDown, BookOpen } from 'lucide-react'

export function Hero() {
  return (
    <section className="relative flex min-h-screen flex-col items-center justify-center overflow-hidden px-6 text-center">
      <div className="absolute inset-0 bg-gradient-to-br from-slate-950 via-slate-900 to-blue-950" />
      <div className="absolute inset-0 bg-[radial-gradient(ellipse_at_top,_rgba(59,130,246,0.15),_transparent_50%)]" />

      <motion.div
        className="relative z-10 max-w-3xl"
        initial={{ opacity: 0, y: 24 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.7 }}
      >
        <motion.p
          className="mb-4 text-sm font-medium uppercase tracking-widest text-blue-400"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.2 }}
        >
          Bioconductor R package
        </motion.p>
        <motion.h1
          className="mb-6 text-5xl font-bold tracking-tight text-white md:text-6xl"
          initial={{ opacity: 0, y: 16 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.3 }}
        >
          evaluomeR
        </motion.h1>
        <motion.p
          className="mb-10 text-lg text-slate-300 md:text-xl"
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          transition={{ delay: 0.45 }}
        >
          Discover the optimal number of clusters for your bioinformatics metrics
          — using bootstrap stability and silhouette quality on real cancer cell-line data.
        </motion.p>

        <motion.div
          className="flex flex-wrap items-center justify-center gap-4"
          initial={{ opacity: 0, y: 12 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.6 }}
        >
          <a
            href="#playground"
            className="inline-flex items-center gap-2 rounded-lg bg-blue-600 px-6 py-3 text-sm font-semibold text-white shadow-lg shadow-blue-600/25 transition hover:bg-blue-500"
          >
            Open Playground
            <ArrowDown className="h-4 w-4" />
          </a>
          <a
            href="https://neobernad.github.io/evaluomeR/"
            className="inline-flex items-center gap-2 rounded-lg border border-slate-700 px-6 py-3 text-sm font-semibold text-slate-200 transition hover:border-slate-500 hover:bg-slate-900"
          >
            <BookOpen className="h-4 w-4" />
            Read Docs
          </a>
        </motion.div>
      </motion.div>

      <motion.p
        className="absolute bottom-8 text-xs text-slate-500"
        initial={{ opacity: 0 }}
        animate={{ opacity: 1 }}
        transition={{ delay: 1 }}
      >
        Precomputed demo · NCI-60 cancer cell lines · k = 3–8 · bootstrap = 100
      </motion.p>
    </section>
  )
}
